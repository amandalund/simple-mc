#include "simple_mc.h"
#include "global.h"

double timer(void)
{
#ifdef _OPENMP
  return omp_get_wtime();
#else
  struct timeval time;
  gettimeofday(&time, 0);
  return time.tv_sec + time.tv_usec/1000000.0;
#endif
}

void copy_particle(Particle *dest, Particle *source)
{
  dest->alive = source->alive;
  dest->energy = source->energy;
  dest->last_energy = source->last_energy;
  dest->mu = source->mu;
  dest->phi = source->phi;
  dest->u = source->u;
  dest->v = source->v;
  dest->w = source->w;
  dest->x = source->x;
  dest->y = source->y;
  dest->z = source->z;
  dest->event = source->event;

  return;
}

