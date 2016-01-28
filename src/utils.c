#include "header.h"

// Linear congruential random number generator
double rn(unsigned long *seed)
{
  unsigned long a = 16807;
  unsigned long m = 2147483647;

  *seed = (a*(*seed)) % m;

  return (double) *seed/m;
}

// Linear congruential random number generator for integer in range [min max)
int rni(unsigned long *seed, int min, int max)
{
  unsigned long a = 16807;
  unsigned long m = 2147483647;

  *seed = (a*(*seed)) % m;

  return min + (int) (max*(*seed)/(m + min));
}

double timer(void)
{
  struct timeval time;
  gettimeofday(&time, 0);
  return time.tv_sec + time.tv_usec/1000000.0;
}
