#include "header.h"

void score_tally(Tally *t, Particle *p)
{
  int ix, iy;

  // Find the indices of the grid box of the particle
  ix = p->x/t->dx;
  iy = p->y/t->dy;

  // Increment number of collisions in this mesh element
  t->sum[ix + t->n*iy]++;

  return;
}

// Simple flux tally
void batch_tally(Tally *t, Parameters *params)
{
  int i, n;
  double xs_t;
  double vol;
  double scale;

  // Number of grid boxes
  n = t->n * t->n;

  // Total cross section of material
  xs_t = params->xs_a + params->xs_f + params->xs_s;

  // Volume
  vol = t->dx * t->dy;

  // For calculating sample estimate of scalar flux
  scale = 1./(vol * xs_t * params->n_particles);

  // Estimate of scalar flux
  for(i=0; i<n; i++){
    t->mean[i] = t->sum[i] * scale;
  }

  // Zero out for next batch
  memset(t->sum, 0, n*sizeof(int));

  return;
}
