#include "header.h"

void score_tally(Tally *t, Particle *p)
{
  int ix, iy, iz;

  // Find the indices of the grid box of the particle
  ix = p->x/t->dx;
  iy = p->y/t->dy;
  iz = p->z/t->dz;

  // Increment number of collisions in this mesh element
  t->sum[ix + t->n*iy + t->n*t->n*iz]++;

  return;
}

// Simple flux tally
void batch_tally(Tally *t, Parameters *params)
{
  int i, n;
  double macro_xs_t;
  double vol;
  double scale;

  // Number of grid boxes
  n = t->n * t->n * t->n;

  // Total cross section of material
  macro_xs_t = params->macro_xs_a + params->macro_xs_f + params->macro_xs_e;

  // Volume
  vol = t->dx * t->dy * t->dz;

  // For calculating sample estimate of scalar flux
  scale = 1./(vol * macro_xs_t * params->n_particles);

  // Estimate of scalar flux
  for(i=0; i<n; i++){
    t->mean[i] = t->sum[i] * scale;
  }

  // Zero out for next batch
  memset(t->sum, 0, n*sizeof(int));

  return;
}
