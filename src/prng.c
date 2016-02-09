#include "header.h"

typedef struct RNG_Parameters_{
  unsigned long long mult;   // multiplier g
  unsigned long long mod;    // modulus, 2^M
  unsigned long long inc;    // increment c
  unsigned long long stride; // stride
  unsigned long long mask;   // mask, 2^M - 1
  unsigned long long period; // period, 2^M for c!=0, else 2^(M-2)
} RNG_Parameters;

// LCG parameters from 'The MCNP5 Random Number Generator', Forrest Brown
// Additional mult for M=63: 2806196910506780709ULL  3249286849523012805ULL
//static const RNG_Parameters RNG = {19073486328125ULL, 281474976710656ULL, 1ULL, 152917, 281474976710655ULL, 70368744177664ULL}; // period 2^46
static const RNG_Parameters RNG = {9219741426499971445ULL, 9223372036854775808ULL, 1ULL, 152917, 9223372036854775807ULL, 9223372036854775808ULL}; // period 2^63

int stream;
unsigned long long seed0[N_STREAMS];
unsigned long long seed[N_STREAMS];

// Linear congruential random number generator: seed = (mult*seed + inc) % mod
double rn(void)
{
  seed[stream] = (RNG.mult*seed[stream] + RNG.inc) & RNG.mask;

  return (double) seed[stream]/RNG.mod;
}

// Linear congruential random number generator for integer in range [a b)
int rni(int a, int b)
{
  seed[stream] = (RNG.mult*seed[stream] + RNG.inc) & RNG.mask;

  return a + (int) (b*seed[stream]/(RNG.mod + a));
}

// Set random number stream
void set_stream(int rn_stream)
{
  stream = rn_stream;

  return;
}

// Set inital seed for each random number sequence
void set_initial_seed(unsigned long long rn_seed0)
{
  int i;

  for(i=0; i<N_STREAMS; i++){
    seed0[i] = rn_seed0 + i;
    seed[i] = seed0[i];
  }

  return;
}

// Algorithm to skip ahead n*RNG.stride random numbers in O(log2(n)) operation,
// from 'The MCNP5 Random Number Generator', Forrest Brown, LA-UR-07K-7961.
void rn_skip(long long n)
{
  unsigned long long g = RNG.mult;
  unsigned long long c = RNG.inc;
  unsigned long long g_new = 1;
  unsigned long long c_new = 0;

  n *= RNG.stride;

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

  seed[stream] = (g_new*seed0[stream] + c_new) & RNG.mask;
}
