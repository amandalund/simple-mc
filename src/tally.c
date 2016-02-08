#include "header.h"
#include "global.h"

// Simple flux tally
void score_tally(Particle *p)
{
  int ix, iy, iz;
  double vol;

  // Volume
  vol = t->dx * t->dy * t->dz;

  // Find the indices of the grid box of the particle
  ix = p->x/t->dx;
  iy = p->y/t->dy;
  iz = p->z/t->dz;

#pragma omp atomic
  // Scalar flux
  t->flux[ix + t->n*iy + t->n*t->n*iz] += 1./(vol * m->xs_t * params->n_particles);

  return;
}
