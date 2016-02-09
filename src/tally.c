#include "header.h"

// Simple flux tally
void score_tally(Parameters *parameters, Material *material, Tally *tally, Particle *p)
{
  int ix, iy, iz;
  double vol;

  // Volume
  vol = tally->dx * tally->dy * tally->dz;

  // Find the indices of the grid box of the particle
  ix = p->x/tally->dx;
  iy = p->y/tally->dy;
  iz = p->z/tally->dz;

  // Scalar flux
  tally->flux[ix + tally->n*iy + tally->n*tally->n*iz] += 1./(vol * material->xs_t * parameters->n_particles);

  return;
}

void reset_tally(Tally *tally)
{
  memset(tally->flux, 0, tally->n*tally->n*tally->n*sizeof(double));

  return;
}
