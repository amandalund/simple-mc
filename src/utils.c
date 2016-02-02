#include "header.h"

// Linear congruential random number generator: seed = (mult*seed + inc) % mod
double rn(unsigned long long *seed)
{
  *seed = (RNG.mult*(*seed) + RNG.inc) & RNG.mask;

  return (double) *seed/RNG.mod;
}

// Linear congruential random number generator for integer in range [a b)
int rni(unsigned long long *seed, int a, int b)
{
  *seed = (RNG.mult*(*seed) + RNG.inc) & RNG.mask;

  return a + (int) (b*(*seed)/(RNG.mod + a));
}

// Algorithm to skip ahead n random numbers in O(log2(n)) operation, from 'The
// MCNP5 Random Number Generator', Forrest Brown, LA-UR-07K-7961.
unsigned long long rn_skip(unsigned long long seed, long long n)
{
  unsigned long long g = RNG.mult;
  unsigned long long c = RNG.inc;
  unsigned long long g_new = 1;
  unsigned long long c_new = 0;

  // Add period until greater than 0
  while(n < 0) n += RNG.period;

  // n % mod
  n = n & RNG.mask;

  // Get mult = mult^n in log2(n) operations
  while(n > 0){
    if(n & 1){
      g_new = g_new*g & RNG.mask;
      c_new = (c_new*g + c) & RNG.mask;
    }
    c = (c*g + c) & RNG.mask;
    g = g*g & RNG.mask;
    n >>= 1;
  }

  return (g_new*seed + c_new) & RNG.mask;
}

double timer(void)
{
  struct timeval time;
  gettimeofday(&time, 0);
  return time.tv_sec + time.tv_usec/1000000.0;
}
