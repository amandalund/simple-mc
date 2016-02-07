#include "header.h"

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
