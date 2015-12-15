#include "header.h"

// Simple flux tally
void score_tally(Tally *t, Particle *p, Material *m, Parameters *params)
{
  int ix, iy;
  double vol;

  // Volume
  vol = t->dx * t->dy;

  // Find the indices of the grid box of the particle
  ix = p->x/t->dx;
  iy = p->y/t->dy;

  // Scalar flux
  t->flux[ix + t->n*iy] += 1./(vol * m->xs_t * params->n_particles);

  return;
}
