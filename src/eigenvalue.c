#include "header.h"

void converge_source(Parameters *params, Bank *source_bank, Bank *fission_bank, Geometry *g, Material *m, Tally *t, Queue *delay_queue)
{
  int i_b;            // index over batches
  int i_g;            // index over generations
  unsigned long i_p;  // index over particles
  double keff_batch;  // keff of batch
  double keff_gen = 1;// keff of generation
  double H;           // shannon entropy
  FILE *fp = NULL;    // file pointer for output

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
        transport(&(source_bank->p[i_p]), g, m, t, fission_bank, keff_gen, params, IGNORE_DB, delay_queue);
      }

      // Calculate generation k_effective and accumulate batch k_effective
      keff_gen = (double) fission_bank->n / source_bank->n;
      keff_batch += keff_gen;

      // Sample new source particles from the particles that were added to the
      // fission bank during this generation
      synchronize_bank(source_bank, fission_bank, g, params);

      // Calculate shannon entropy to assess source convergence
      H = shannon_entropy(g, source_bank, params);
      if(params->write_entropy == TRUE){
        write_entropy(H, fp, params->entropy_file);
      }

      // Write the source distribution
      if(params->write_source == TRUE){
        write_source(g, source_bank, params, fp, params->source_file);
      }
    }

    // Calculate k_effective
    keff_batch /= params->n_generations;

    // Status text
    printf("%-15d %-15f %-15f\n", i_b+1, H, keff_batch);
  }

  return;
}

void build_delay_bank(Parameters *params, Bank *source_bank, Bank *fission_bank, Geometry *g, Material *m, Tally *t, Queue *delay_queue)
{
  int i_g;            // index over generations
  unsigned long i_p;  // index over particles
  double keff_gen = 1;// keff of generation
  double H;           // shannon entropy
  FILE *fp = NULL;    // file pointer for output

  // Loop over generations
  for(i_g=0; i_g<params->lag; i_g++){

    // Loop over particles
    for(i_p=0; i_p<source_bank->n; i_p++){

      // Transport the next particle from source bank
      transport(&(source_bank->p[i_p]), g, m, t, fission_bank, keff_gen, params, BUILD_DB, delay_queue);
    }

    // Calculate generation k_effective and accumulate batch k_effective
    keff_gen = (double) fission_bank->n / source_bank->n;

    // Sample new source particles from the particles that were added to the
    // fission bank during this generation
    synchronize_bank(source_bank, fission_bank, g, params);

    // Calculate shannon entropy to assess source convergence
    H = shannon_entropy(g, source_bank, params);
    if(params->write_entropy == TRUE){
      write_entropy(H, fp, params->entropy_file);
    }
  }

  return;
}

void run_eigenvalue(Parameters *params, Bank *source_bank, Bank *fission_bank, Geometry *g, Material *m, Tally *t, double *keff, Queue *delay_queue)
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
        transport(&(source_bank->p[i_p]), g, m, t, fission_bank, keff_gen, params, USE_DB, delay_queue);
      }

      // Calculate generation k_effective and accumulate batch k_effective
      keff_gen = (double) fission_bank->n / source_bank->n;
      keff_batch += keff_gen;

      // Sample new source particles from the particles that were added to the
      // fission bank during this generation
      synchronize_bank(source_bank, fission_bank, g, params);

      // Calculate shannon entropy to assess source convergence
      H = shannon_entropy(g, source_bank, params);
      if(params->write_entropy == TRUE){
        write_entropy(H, fp, params->entropy_file);
      }

      // Write the source distribution
      if(params->write_source == TRUE){
        write_source(g, source_bank, params, fp, params->source_file);
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

void synchronize_bank(Bank *source_bank, Bank *fission_bank, Geometry *g, Parameters *params)
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
      j = rni(&(params->seed), 0, i+1);
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
      j = rni(&(params->seed), 0, n_f);
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
