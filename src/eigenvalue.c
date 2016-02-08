#include "header.h"
#include "global.h"

void run_eigenvalue(void)
{
  int i_b; // index over batches
  int i_a = -1; // index over active batches
  int i_g; // index over generations
  unsigned long i_p; // index over particles
  double keff_gen = 1; // keff of generation
  double keff_batch; // keff of batch
  double keff_mean; // keff mean over active batches
  double keff_std; // keff standard deviation over active batches
  double H; // shannon entropy
  Particle p;

  // Loop over batches
  for(i_b=0; i_b<params->n_batches; i_b++){

    keff_batch = 0;

    // Write coordinates of particles in source bank
    if(params->write_bank == TRUE){
      write_bank(source_bank, params->bank_file);
    }

    // Turn on tallying and increment index in active batches
    if(i_b >= params->n_batches - params->n_active){
      i_a++;
      if(params->tally == TRUE){
        t->tallies_on = TRUE;
      }
    }

    // Loop over generations
    for(i_g=0; i_g<params->n_generations; i_g++){

      // Set RNG stream for tracking
      set_stream(STREAM_TRACK);

#ifdef _OPENMP
#pragma omp parallel default(none) private(i_p, p) shared(i_b, i_g, params, source_bank, g, m, t)
{
#pragma omp for
#endif
      // Loop over particles
      for(i_p=0; i_p<params->n_particles; i_p++){

	// Set seed for particle i_p by skipping ahead in the random number
	// sequence stride*(total particles simulated) numbers. This allows for
	// reproducibility of the particle history.
        rn_skip((i_b*params->n_generations + i_g)*params->n_particles + i_p);

        // Copy next particle into p
        copy_particle(&p, &(source_bank->p[i_p]));

        // Transport the next particle from source bank
	transport(&p);
      }
#ifdef _OPENMP
}
#endif
      set_stream(STREAM_OTHER);
      rn_skip(i_b*params->n_generations + i_g);

#ifdef _OPENMP
      // Merge fission banks from threads
      merge_fission_banks();
#endif

      // Calculate generation k_effective and accumulate batch k_effective
      keff_gen = (double) fission_bank->n / source_bank->n;
      keff_batch += keff_gen;

      // Sample new source particles from the particles that were added to the
      // fission bank during this generation
      synchronize_bank();

      // Calculate shannon entropy to assess source convergence
      H = shannon_entropy(source_bank);
      if(params->write_entropy == TRUE){
        write_entropy(H, params->entropy_file);
      }

      // Write the source distribution
      if(params->write_source == TRUE){
        write_source(source_bank, params->source_file);
      }
    }

    // Calculate k_effective
    keff_batch /= params->n_generations;
    if(i_a >= 0){
      keff[i_a] = keff_batch;
    }

    // Tallies for this realization
    if(t->tallies_on == TRUE){
      if(params->write_tally == TRUE){
        write_tally(params->tally_file);
      }
      memset(t->flux, 0, t->n*t->n*t->n*sizeof(double));
    }

    // Calculate keff mean and standard deviation
    calculate_keff(&keff_mean, &keff_std, i_a+1);

    // Status text
    if(i_a < 0){
      printf("%-15d %-15f %-15f\n", i_b+1, H, keff_batch);
    }
    else{
      printf("%-15d %-15f %-15f %f +/- %-15f\n", i_b+1, H, keff_batch, keff_mean, keff_std);
    }
  }

  // Write out keff
  if(params->write_keff == TRUE){
    write_keff(keff, params->n_active, params->keff_file);
  }

  if(params->save_source == TRUE){
    save_source(source_bank);
  }

  return;
}

void merge_fission_banks(void)
{
#ifdef _OPENMP
  unsigned long n_sites = 0; // total number of source sites in fission bank
  int i_t; // index over threads

#pragma omp parallel
{
#pragma omp for ordered
  for(i_t=0; i_t<n_threads; i_t++){
#pragma omp ordered
{
    memcpy(&(master_fission_bank->p[n_sites]), fission_bank->p, fission_bank->n*sizeof(Particle));
    n_sites += fission_bank->n;
}
  }
#pragma omp barrier
  // Copy shared fission bank sites into master thread's fission bank
  if(thread_id == 0){
    memcpy(fission_bank->p, master_fission_bank->p, n_sites*sizeof(Particle));
    fission_bank->n = n_sites;
  }
  else{
    fission_bank->n = 0;
  }
}
#endif

  return;
}

void synchronize_bank(void)
{
  unsigned long i, j;
  unsigned long n_s = source_bank->n;
  unsigned long n_f = fission_bank->n;

  // If the fission bank is larger than the source bank, randomly select
  // n_particles sites from the fission bank to create the new source bank
  if(n_f >= n_s){

    // Copy first n_particles sites from fission bank to source bank
    memcpy(source_bank->p, fission_bank->p, n_s*sizeof(Particle));

    // Replace elements with decreasing probability, such that after final
    // iteration each particle in fission bank will have equal probability of
    // being selected for source bank
    for(i=n_s; i<n_f; i++){
      j = rni(0, i+1);
      if(j<n_s){
        memcpy(&(source_bank->p[j]), &(fission_bank->p[i]), sizeof(Particle));
      }
    }
  }

  // If the fission bank is smaller than the source bank, use all fission bank
  // sites for the source bank and randomly sample remaining particles from
  // fission bank
  else{

    // First randomly sample particles from fission bank
    for(i=0; i<(n_s-n_f); i++){
      j = rni(0, n_f);
      memcpy(&(source_bank->p[i]), &(fission_bank->p[j]), sizeof(Particle));
    }

    // Fill remaining source bank sites with fission bank
    memcpy(&(source_bank->p[n_s-n_f]), fission_bank->p, n_f*sizeof(Particle));
  }

  fission_bank->n = 0;

  return;
}

// Calculates the shannon entropy of the source distribution to assess
// convergence
double shannon_entropy(Bank *b)
{
  unsigned long i;
  double H = 0.0;
  double dx, dy, dz;
  unsigned long ix, iy, iz;
  unsigned long n;
  unsigned long *count;
  Particle *p;

  // Determine an appropriate number of grid boxes in each dimension
  n = ceil(pow(b->n/20, 1.0/3.0));

  // Find grid spacing
  dx = g->x/n;
  dy = g->y/n;
  dz = g->z/n;

  // Allocate array to keep track of number of sites in each grid box
  count = calloc(n*n*n, sizeof(unsigned long));

  for(i=0; i<b->n; i++){
    p = &(b->p[i]);

    // Find the indices of the grid box of the particle
    ix = p->x/dx;
    iy = p->y/dy;
    iz = p->z/dz;

    count[ix*n*n + iy*n + iz]++;
  }

  // Calculate the shannon entropy
  for(i=0; i<n*n*n; i++){
    if(count[i] > 0){
      H -= ((double)count[i]/b->n) * log2((double)count[i]/b->n);
    }
  }

  free(count);

  return H;
}

void calculate_keff(double *mean, double *std, int n)
{
  int i;

  *mean = 0;
  *std = 0;

  // Calculate mean
  for(i=0; i<n; i++){
    *mean += keff[i];
  }
  *mean /= n;

  // Calculate standard deviation
  for(i=0; i<n; i++){
    *std += pow(keff[i] - *mean, 2);
  }
  *std = sqrt(*std/(n-1));

  return;
}
