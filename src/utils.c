#include "simple_mc.h"

double timer(void)
{
  struct timeval time;
  gettimeofday(&time, 0);
  return time.tv_sec + time.tv_usec/1000000.0;
}

void copy_particle(Particle *dest, Particle *source)
{
  dest->alive = source->alive;
  dest->mu = source->mu;
  dest->phi = source->phi;
  dest->u = source->u;
  dest->v = source->v;
  dest->w = source->w;
  dest->x = source->x;
  dest->y = source->y;
  dest->z = source->z;

  return;
}

