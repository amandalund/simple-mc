#include "simple_mc.h"

void ramp_up(Parameters *parameters, Geometry *geometry, Material *material, Bank *source_bank, Bank *fission_bank, Tally *tally)
{
  int i_s; // index over stages
  int i_g; // index over generations
  unsigned long i_p;  // index over particles
  double keff_gen = 1; // keff of generation
  double keff_stage; // keff of stage
  double H; // shannon entropy
  Particle p;

  printf("%-15s %-15s %-15s\n", "STAGE", "ENTROPY", "KEFF");

  // Loop over stages
  for(i_s=0; i_s<parameters->cnvg_n_stages; i_s++){

    // Switch RNG stream off tracking
    set_stream(STREAM_OTHER);

    keff_stage = 0;

    add_particles(parameters, geometry, source_bank, i_s);

    // Loop over generations
    for(i_g=0; i_g<parameters->cnvg_n_generations[i_s]; i_g++){

      // Set RNG stream for tracking
      set_stream(STREAM_TRACK);

      // Loop over particles
      for(i_p=0; i_p<source_bank->n; i_p++){

	// Set seed for particle i_p by skipping ahead in the random number
	// sequence stride*(total particles simulated) numbers from the initial
	// seed. This allows for reproducibility of the particle history.
        rn_skip((i_s*parameters->n_generations + i_g)*parameters->n_particles + i_p);

        // Copy next particle into p
        copy_particle(&p, &(source_bank->p[i_p]));

        // Transport the next particle
        transport(parameters, geometry, material, source_bank, fission_bank, tally, &p);

        // Increment the total number of histories simulated
        parameters->n_histories++;
      }

      // Switch RNG stream off tracking
      set_stream(STREAM_OTHER);

      // Calculate generation k_effective and accumulate batch k_effective
      keff_gen = (double) fission_bank->n / source_bank->n;
      keff_stage += keff_gen;

      // Sample new source particles from the particles that were added to the
      // fission bank during this generation
      synchronize_bank(source_bank, fission_bank);

      // Calculate shannon entropy to assess source convergence
      H = shannon_entropy(geometry, source_bank, parameters);
      if(parameters->write_entropy == TRUE){
        write_entropy(H, parameters->entropy_file);
      }

      // Write the histories
      if(parameters->write_histories == TRUE){
        write_histories(parameters->n_histories, parameters->histories_file);
      }
    }

    // Calculate k_effective
    keff_stage /= parameters->cnvg_n_generations[i_s];

    printf("%-15d %-15f %-15f\n", i_s+1, H, keff_stage);
  }

  printf("Converged after %lu histories.\n", parameters->n_histories);

  return;
}

void add_particles(Parameters *parameters, Geometry *geometry, Bank *source_bank, int i_s)
{
  unsigned long i_p;  // index over particles
  int n; // total number of particles in stage
  int n_prev; // total number of particles in previousstage
  int i;

  // the total number of particles in the source bank at this stage
  n = parameters->cnvg_n_particles[i_s];

  // at first stage sample initial source particles
  if(i_s == 0){
    for(i_p=0; i_p<n; i_p++){
      sample_source_particle(geometry, &(source_bank->p[i_p]));
      source_bank->n++;
    }
  }

  // At remaining stages add new particles
  else{

    n_prev = parameters->cnvg_n_particles[i_s-1];

    // If the number of particles in the previous stage is greater than the
    // number to be added in this stage, randomly select n - n_prev sites from
    // the previous stage
    if(n_prev >= n - n_prev){

      // Copy first n - n_prev sites from previous stage
      memcpy(&(source_bank->p[n_prev]), source_bank->p, (n-n_prev)*sizeof(Particle));

      // Replace elements with decreasing probability, such that after final
      // iteration each particle from previous stage will have equal probability
      // of being selected
      for(i_p=n-n_prev; i_p<n_prev; i_p++){
        i = rni(0, i_p+1);
        if(i<n-n_prev){
          memcpy(&(source_bank->p[n_prev+i]), &(source_bank->p[i_p]), sizeof(Particle));
        }
      }
    }

    // If the number of particles in the previous stage is less than the number
    // to be added to this stage, use all particles from previous stage and
    // randomly sample remaining particles from previous stage
    else{

      // First copy particles from previous stage
      memcpy(&(source_bank->p[n_prev]), source_bank->p, n_prev*sizeof(Particle));

      // Fill in remaining particles by sampling from previous stage
      for(i_p=2*n_prev; i_p<n; i_p++){
        i = rni(0, n_prev);
        memcpy(&(source_bank->p[i_p]), &(source_bank->p[i]), sizeof(Particle));
      }
    }

    source_bank->n = n;
  }

  return;
}

void run_eigenvalue(Parameters *parameters, Geometry *geometry, Material *material, Bank *source_bank, Bank *fission_bank, Tally *tally, double *keff)
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

  printf("%-15s %-15s %-15s %-15s\n", "BATCH", "ENTROPY", "KEFF", "MEAN KEFF");

  // Loop over batches
  for(i_b=0; i_b<parameters->n_batches; i_b++){

    keff_batch = 0;

    // Write coordinates of particles in source bank
    if(parameters->write_bank == TRUE){
      write_bank(source_bank, parameters->bank_file);
    }

    // Turn on tallying and increment index in active batches
    if(i_b >= parameters->n_batches - parameters->n_active){
      i_a++;
      if(parameters->tally == TRUE){
        tally->tallies_on = TRUE;
      }
    }

    // Loop over generations
    for(i_g=0; i_g<parameters->n_generations; i_g++){

      // Set RNG stream for tracking
      set_stream(STREAM_TRACK);

      // Loop over particles
      for(i_p=0; i_p<parameters->n_particles; i_p++){

	// Set seed for particle i_p by skipping ahead in the random number
	// sequence stride*(total particles simulated) numbers from the initial
	// seed. This allows for reproducibility of the particle history.
        rn_skip((i_b*parameters->n_generations + i_g)*parameters->n_particles + i_p);

        // Copy next particle into p
        copy_particle(&p, &(source_bank->p[i_p]));

        // Transport the next particle
        transport(parameters, geometry, material, source_bank, fission_bank, tally, &p);

        // Increment the total number of histories simulated
        parameters->n_histories++;
      }

      // Switch RNG stream off tracking
      set_stream(STREAM_OTHER);

      // Calculate generation k_effective and accumulate batch k_effective
      keff_gen = (double) fission_bank->n / source_bank->n;
      keff_batch += keff_gen;

      // Sample new source particles from the particles that were added to the
      // fission bank during this generation
      synchronize_bank(source_bank, fission_bank);

      // Calculate shannon entropy to assess source convergence
      H = shannon_entropy(geometry, source_bank, parameters);
      if(parameters->write_entropy == TRUE){
        write_entropy(H, parameters->entropy_file);
      }

      // Write the total number of histories run so far
      if(parameters->write_histories == TRUE){
        write_histories(parameters->n_histories, parameters->histories_file);
      }

      // Write the source distribution
      if(parameters->write_source == TRUE){
        write_source(parameters, geometry, source_bank, parameters->source_file);
      }
    }

    // Calculate k_effective
    keff_batch /= parameters->n_generations;
    if(i_a >= 0){
      keff[i_a] = keff_batch;
    }

    // Tallies for this realization
    if(tally->tallies_on == TRUE){
      if(parameters->write_tally == TRUE){
        write_tally(tally, parameters->tally_file);
      }
      memset(tally->flux, 0, tally->n*tally->n*tally->n*sizeof(double));
    }

    // Calculate keff mean and standard deviation
    calculate_keff(keff, &keff_mean, &keff_std, i_a+1);

    // Status text
    if(i_a < 0){
      printf("%-15d %-15f %-15f\n", i_b+1, H, keff_batch);
    }
    else{
      printf("%-15d %-15f %-15f %f +/- %-15f\n", i_b+1, H, keff_batch, keff_mean, keff_std);
    }
  }

  // Write out keff
  if(parameters->write_keff == TRUE){
    write_keff(keff, parameters->n_active, parameters->keff_file);
  }

  if(parameters->save_source == TRUE){
    save_source(source_bank);
  }

  return;
}

void synchronize_bank(Bank *source_bank, Bank *fission_bank)
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
double shannon_entropy(Geometry *geometry, Bank *b, Parameters *parameters)
{
  unsigned long i;
  double H = 0.0;
  double dx, dy, dz;
  unsigned long ix, iy, iz;
  unsigned long n;
  unsigned long *count;
  Particle *p;

  // Determine an appropriate number of grid boxes in each dimension
  n = parameters->entropy_bins;

  // Find grid spacing
  dx = geometry->Lx/n;
  dy = geometry->Ly/n;
  dz = geometry->Lz/n;

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

void calculate_keff(double *keff, double *mean, double *std, int n)
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
