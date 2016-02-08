#include "header.h"
#include "global.h"

// Simple flux tally
void score_tally(Tally *t, Particle *p)
{
  int ix, iy, iz;
  double vol;

  // Volume
  vol = t->dx * t->dy * t->dz;

  // Find the indices of the grid box of the particle
  ix = p->x/t->dx;
  iy = p->y/t->dy;
  iz = p->z/t->dz;

#ifdef _OPENMP
#pragma omp atomic
#endif
  // Scalar flux
  t->flux[ix + t->n*iy + t->n*t->n*iz] += 1./(vol * material->xs_t * parameters->n_particles);

  return;
}
