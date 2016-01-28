#include "header.h"

// Linear congruential random number generator
double rn(unsigned long long *seed)
{
  unsigned long long a = 19073486328125;
  unsigned long long m = 281474976710656;

  *seed = (a*(*seed)) % m;

  return (double) *seed/m;
}

// Linear congruential random number generator for integer in range [min max)
int rni(unsigned long long *seed, int min, int max)
{
  unsigned long long a = 19073486328125;
  unsigned long long m = 281474976710656;

  *seed = (a*(*seed)) % m;

  return min + (int) (max*(*seed)/(m + min));
}

double timer(void)
{
  struct timeval time;
  gettimeofday(&time, 0);
  return time.tv_sec + time.tv_usec/1000000.0;
}
