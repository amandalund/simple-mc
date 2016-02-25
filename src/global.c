#include "simple_mc.h"
#include "global.h"

Bank *fission_bank;
#ifdef _OPENMP
Bank *master_fission_bank;
int n_threads;
int thread_id;
#pragma omp threadprivate(fission_bank, thread_id)
#endif
