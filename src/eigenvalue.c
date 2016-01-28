#include "header.h"

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
