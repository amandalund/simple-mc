#include "header.h"

void ramp_up(Parameters *params, Bank *source_bank, Bank *fission_bank, Geometry *g, Material *m, Tally *t)
{
  int i_s;            // index over stages
  int i_g;            // index over generations
  unsigned long i_p;  // index over particles
  double keff_gen = 1;// keff of generation
  double H;           // shannon entropy
  FILE *fp = NULL;    // file pointer for output
  int n;              // total number of particles in stage
  int n_prev;         // total number of particles in previousstage
  int i;

  // Loop over stages
  for(i_s=0; i_s<params->cnvg_n_stages; i_s++){

    n = params->cnvg_n_particles[i_s];

    // At first stage sample initial source particles
    if(i_s == 0){
      for(i_p=0; i_p<n; i_p++){
        sample_source_particle(&(source_bank->p[i_p]), g);
        source_bank->n++;
      }
    }

    // At remaining stages add new particles
    else{

      n_prev = params->cnvg_n_particles[i_s-1];

      // If the number of particles in the previous stage is greater than the
      // number to be added in this stage, randomly select n - n_prev sites from
      // the previous stage
      if(n_prev >= n - n_prev){

        // Copy first n - n_prev sites from previous stage
        memcpy(&(source_bank->p[n_prev]), source_bank->p, (n - n_prev)*sizeof(Particle));

        // Replace elements with decreasing probability, such that after final
        // iteration each particle from previous stage will have equal probability
        // of being selected
        for(i_p=n-n_prev; i_p<n_prev; i_p++){
          i = rand() % (i_p+1);
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
          i = rand() % n_prev;
          memcpy(&(source_bank->p[i_p]), &(source_bank->p[i]), sizeof(Particle));
        }
      }

      source_bank->n = n;
    }

    // Write coordinates of particles in source bank
    if(params->write_bank == TRUE){
      write_bank(source_bank, fp, params->bank_file);
    }

    // Loop over generations
    for(i_g=0; i_g<params->cnvg_n_generations[i_s]; i_g++){

      // Loop over particles
      for(i_p=0; i_p<source_bank->n; i_p++){

        // Transport the next particle from source bank
        transport(&(source_bank->p[i_p]), g, m, t, fission_bank, keff_gen, params);

        // Increment total number of histories
//        params->n_histories++;
      }

      // Calculate generation k_effective and accumulate batch k_effective
      keff_gen = (double) fission_bank->n / source_bank->n;

      // Sample new source particles from the particles that were added to the
      // fission bank during this generation
      synchronize_bank(source_bank, fission_bank, g);

      // Calculate shannon entropy to assess source convergence
      H = shannon_entropy(g, source_bank, params);
      if(params->write_entropy == TRUE){
        write_entropy(H, fp, params->entropy_file);
      }
    }
  }

  return;
}

/*void ramp_up(Parameters *params, Bank *source_bank, Bank *fission_bank, Geometry *g, Material *m, Tally *t)
{
  int i_s;            // index over stages
  int i_g;            // index over generations
  unsigned long i_p;  // index over particles
  double keff_gen = 1;// keff of generation
  double H;           // shannon entropy
  FILE *fp = NULL;    // file pointer for output
  int n_particles;    // total number of particles in stage

  // Loop over stages
  for(i_s=0; i_s<params->cnvg_n_stages; i_s++){

    n_particles = params->cnvg_n_particles[i_s];

    // At first stage sample initial source particles
    if(i_s == 0){
      for(i_p=0; i_p<n_particles; i_p++){
        sample_source_particle(&(source_bank->p[i_p]), g);
        source_bank->n++;
      }
    }

    // At remaining stages add new uniformly sampled particles
    else{

      // Iterate over particles to add and  uniformly sample
      for(i_p=params->cnvg_n_particles[i_s-1]; i_p<n_particles; i_p++){

        // Sample particle uniformly from within bin
        sample_source_particle(&(source_bank->p[i_p]), g);
        source_bank->n++;
      }
    }

    // Write coordinates of particles in source bank
    if(params->write_bank == TRUE){
      write_bank(source_bank, fp, params->bank_file);
    }

    // Loop over generations
    for(i_g=0; i_g<params->cnvg_n_generations[i_s]; i_g++){

      // Loop over particles
      for(i_p=0; i_p<source_bank->n; i_p++){

        // Transport the next particle from source bank
        transport(&(source_bank->p[i_p]), g, m, t, fission_bank, keff_gen, params);

        // Increment total number of histories
        params->n_histories++;
      }

      // Calculate generation k_effective and accumulate batch k_effective
      keff_gen = (double) fission_bank->n / source_bank->n;

      // Sample new source particles from the particles that were added to the
      // fission bank during this generation
      synchronize_bank(source_bank, fission_bank, g);

      // Calculate shannon entropy to assess source convergence
      H = shannon_entropy(g, source_bank, params);
      if(params->write_entropy == TRUE){
        write_entropy(H, fp, params->entropy_file);
      }
    }
  }

  return;
}

void ramp_up(Parameters *params, Bank *source_bank, Bank *fission_bank, Geometry *g, Material *m, Tally *t)
{
  int i_s;            // index over stages
  int i_g;            // index over generations
  int i_m;            // index over mesh
  unsigned long i_p;  // index over particles
  double keff_gen = 1;// keff of generation
  double H;           // shannon entropy
  FILE *fp = NULL;    // file pointer for output
  int n_bins;         // number of bins in each dimension
  int n;              // total number of bins
  int n_particles;    // total number of particles in stage
  int *count;         // fraction of particles in each bin
  double dx, dy, dz;
  int ix, iy, iz;
  Particle *p;
  double cutoff;
  double prob = 0.0;
  int i;

  n_bins = params->cnvg_n_bins;
  n = n_bins*n_bins*n_bins;
  count = malloc(n*sizeof(int));

  // Find grid spacing
  dx = g->x/n_bins;
  dy = g->y/n_bins;
  dz = g->z/n_bins;

  // Loop over stages
  for(i_s=0; i_s<params->cnvg_n_stages; i_s++){

    n_particles = params->cnvg_n_particles[i_s];

    // At first stage sample initial source particles
    if(i_s == 0){
      for(i_p=0; i_p<n_particles; i_p++){
        sample_source_particle(&(source_bank->p[i_p]), g);
        source_bank->n++;
      }
    }

    // At remaining stages add new particles -- count number of particles in
    // each bin, add new uniformly sampled particles to each bin proportional
    // to number of particles that bin contained
    else{

      memset(count, 0, n*sizeof(int));

      // Count number of particles in each bin
      for(i_p=0; i_p<params->cnvg_n_particles[i_s-1]; i_p++){
        p = &(source_bank->p[i_p]);

        // Find the indices of the grid box of the particle
        ix = p->x/dx;
        iy = p->y/dy;
        iz = p->z/dz;

        count[ix*n_bins*n_bins + iy*n_bins + iz]++;
      }

      // Iterate over particles to add, sample bin, then uniformly sample from
      // within that bin
      for(i_p=params->cnvg_n_particles[i_s-1]; i_p<n_particles; i_p++){

        i = 0;

        // Cutoff for sampling bin
        cutoff = rn()*params->cnvg_n_particles[i_s-1];

        // Sample which bin particle goes in
        while(prob < cutoff){
          i_m = i;
          prob += count[i];
          i++;
        }

        iz = i_m % n_bins;
        iy = (i_m/n_bins) % n_bins;
        ix = i_m / (n_bins*n_bins);

        // Sample particle uniformly from within bin
        sample_bounded_source_particle(&(source_bank->p[i_p]), ix*dx, (ix+1)*dx, iy*dy, (iy+1)*dy, iz*dz, (iz+1)*dz);
      }
      source_bank->n = n_particles;
    }

    // Write coordinates of particles in source bank
    if(params->write_bank == TRUE){
      write_bank(source_bank, fp, params->bank_file);
    }

    // Loop over generations
    for(i_g=0; i_g<params->cnvg_n_generations[i_s]; i_g++){

      // Loop over particles
      for(i_p=0; i_p<source_bank->n; i_p++){

        // Transport the next particle from source bank
        transport(&(source_bank->p[i_p]), g, m, t, fission_bank, keff_gen, params);

        // Increment total number of histories
        params->n_histories++;
      }

      // Calculate generation k_effective and accumulate batch k_effective
      keff_gen = (double) fission_bank->n / source_bank->n;

      // Sample new source particles from the particles that were added to the
      // fission bank during this generation
      synchronize_bank(source_bank, fission_bank, g);

      // Calculate shannon entropy to assess source convergence
      H = shannon_entropy(g, source_bank, params);
      if(params->write_entropy == TRUE){
        write_entropy(H, fp, params->entropy_file);
      }
    }
  }

  free(count);

  return;
}
*/
void converge_source(Parameters *params, Bank *source_bank, Bank *fission_bank, Geometry *g, Material *m, Tally *t)
{
  int i_b;            // index over batches
  int i_g;            // index over generations
  unsigned long i_p;  // index over particles
  double keff_batch;  // keff of batch
  double keff_gen = 1;// keff of generation
  double H;           // shannon entropy
  FILE *fp = NULL;    // file pointer for output

  // Sample source particles or load a source
  if(params->load_source == TRUE){
    load_source(source_bank);
    source_bank->n = params->n_particles;
  }
  else{
    for(i_p=0; i_p<params->n_particles; i_p++){
      sample_source_particle(&(source_bank->p[i_p]), g);
      source_bank->n++;
    }
  }

  // Loop over inactive batches
  for(i_b=0; i_b<params->n_batches - params->n_active; i_b++){

    keff_batch = 0;

    // Write coordinates of particles in source bank
    if(params->write_bank == TRUE){
      write_bank(source_bank, fp, params->bank_file);
    }

    // Loop over generations
    for(i_g=0; i_g<params->n_generations; i_g++){

      // Loop over particles
      for(i_p=0; i_p<source_bank->n; i_p++){

        // Transport the next particle from source bank
        transport(&(source_bank->p[i_p]), g, m, t, fission_bank, keff_gen, params);

        // Increment total number of histories
//        params->n_histories++;
      }

      // Calculate generation k_effective and accumulate batch k_effective
      keff_gen = (double) fission_bank->n / source_bank->n;
      keff_batch += keff_gen;

      // Sample new source particles from the particles that were added to the
      // fission bank during this generation
      synchronize_bank(source_bank, fission_bank, g);

      // Calculate shannon entropy to assess source convergence
      H = shannon_entropy(g, source_bank, params);
      if(params->write_entropy == TRUE){
        write_entropy(H, fp, params->entropy_file);
      }
    }

    // Calculate k_effective
    keff_batch /= params->n_generations;

    // Status text
    printf("%-15d %-15f %-15f\n", i_b+1, H, keff_batch);
  }

  return;
}

void run_eigenvalue(Parameters *params, Bank *source_bank, Bank *fission_bank, Geometry *g, Material *m, Tally *t, double *keff)
{
  int i_b;            // index over batches
  int i_g;            // index over generations
  unsigned long i_p;  // index over particles
  FILE *fp = NULL;    // file pointer for output
  double keff_gen = 1;// keff of generation
  double keff_batch;  // keff of batch
  double keff_mean;   // keff mean over active batches
  double keff_std;    // keff standard deviation over active batches
  double H;           // shannon entropy

  // Turn on tallying
  if(params->tally == TRUE){
    t->tallies_on = TRUE;
  }

  // Loop over active batches
  for(i_b=0; i_b<params->n_active; i_b++){

    keff_batch = 0;

    // Write coordinates of particles in source bank
    if(params->write_bank == TRUE){
      write_bank(source_bank, fp, params->bank_file);
    }

    // Loop over generations
    for(i_g=0; i_g<params->n_generations; i_g++){

      // Loop over particles
      for(i_p=0; i_p<source_bank->n; i_p++){

        // Transport the next particle from source bank
        transport(&(source_bank->p[i_p]), g, m, t, fission_bank, keff_gen, params);

        // Increment total number of histories
//        params->n_histories++;
      }

      // Calculate generation k_effective and accumulate batch k_effective
      keff_gen = (double) fission_bank->n / source_bank->n;
      keff_batch += keff_gen;

      // Sample new source particles from the particles that were added to the
      // fission bank during this generation
      synchronize_bank(source_bank, fission_bank, g);

      // Calculate shannon entropy to assess source convergence
      H = shannon_entropy(g, source_bank, params);
      if(params->write_entropy == TRUE){
        write_entropy(H, fp, params->entropy_file);
      }
    }

    // Calculate k_effective
    keff_batch /= params->n_generations;
    keff[i_b] = keff_batch;

    // Tallies for this realization
    if(t->tallies_on == TRUE){
      if(params->write_tally == TRUE){
        write_tally(t, fp, params->tally_file);
      }
      memset(t->flux, 0, t->n*t->n*t->n*sizeof(double));
    }

    // Calculate keff mean and standard deviation
    calculate_keff(keff, &keff_mean, &keff_std, i_b+1);

    // Status text
    printf("%-15d %-15f %-15f %f +/- %-15f\n", params->n_batches - params->n_active + i_b + 1, H, keff_batch, keff_mean, keff_std);
  }

  return;
}

void synchronize_bank(Bank *source_bank, Bank *fission_bank, Geometry *g)
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
      j = rand() % (i+1);
      if(j<n_s){
        memcpy(&(source_bank->p[j]), &(fission_bank->p[i]), sizeof(Particle));
      }
    }
  }

  // If the fission bank is smaller than the source bank, use all fission bank
  // sites for the source bank and randomly sample remaining particles from
  // source distribution
  else{

    // First sample particles from source distribution
    for(i=0; i<(n_s-n_f); i++){
      sample_source_particle(&(source_bank->p[i]), g);
    }

    // Fill remaining source bank sites with fission bank
    memcpy(&(source_bank->p[n_s-n_f]), fission_bank->p, n_f*sizeof(Particle));
  }

  fission_bank->n = 0;

  return;
}

// Calculates the shannon entropy of the source distribution to assess
// convergence
double shannon_entropy(Geometry *g, Bank *b, Parameters *params)
{
  unsigned long i;
  double H = 0.0;
  double dx, dy, dz;
  unsigned long ix, iy, iz;
  unsigned long n;
  unsigned long *count;
  Particle *p;

  // Determine an appropriate number of grid boxes in each dimension
  n = params->entropy_bins;

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
