#include "header.h"
#include "global.h"

Parameters *parameters;
Geometry *geometry;
Material *material;
Tally *tally;
Bank *source_bank;
Bank *fission_bank;
double *keff;

#ifdef _OPENMP
Bank *master_fission_bank;
int n_threads;
int thread_id;
#pragma omp threadprivate(fission_bank, thread_id)
#endif
