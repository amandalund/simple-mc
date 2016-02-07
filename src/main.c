#include "header.h"

Bank *fission_bank;
#ifdef _OPENMP
int thread_id;
#pragma omp threadprivate(fission_bank, thread_id)
#endif

int main(int argc, char *argv[])
{
  unsigned long long seed1, seed2; // RNG seeds for different RN sequences
  unsigned long long seed1_0, seed2_0; // initial RNG seeds
  int i_b;            // index over batches
  int i_a = -1;       // index over active batches
  int i_g;            // index over generations
  unsigned long i_p;  // index over particles
  FILE *fp = NULL;    // file pointer for output
  double t1, t2;      // timers
  double *keff;       // effective multiplication factor
  double keff_gen = 1;// keff of generation
  double keff_batch;  // keff of batch
  double keff_mean;   // keff mean over active batches
  double keff_std;    // keff standard deviation over active batches
  double H;           // shannon entropy
  Particle p;
  Parameters *params; // user defined parameters
  Geometry *g;
  Material *m;
  Tally *t;
  Bank *source_bank;
#ifdef _OPENMP
  int i_t;            // index over threads
  int n_threads;      // number of OpenMP threads
  unsigned long n_sites; // total number of source sites in fission bank
  Bank *master_fission_bank;
#endif

  // Get inputs: set parameters to default values, then parse parameter file,
  // then override with any command line inputs
  params = init_params();
  parse_params(params);
  read_CLI(argc, argv, params);
  print_params(params);

  // Set RNG seeds for two different RN sequences: one for tracking and one for
  // everything else
  seed1 = params->seed;
  seed2 = seed1 + 1;

  // Set up output files
  init_output(params, fp);

  // Set up array for keff
  keff = calloc(params->n_active, sizeof(double));

  // Set up geometry
  g = init_geometry(params);

  // Set up material
  m = init_material(params, &seed2);

  // Set up tallies
  t = init_tally(params);

#ifdef _OPENMP
  // Set number of openmp threads
  omp_set_num_threads(params->n_threads);

  // Allocate one fission bank for each thread and one master fission bank to
  // collect fission sites at end of each generation
#pragma omp parallel
{
  n_threads = omp_get_num_threads();
  thread_id = omp_get_thread_num();
  if(thread_id == 0){
    fission_bank = init_bank(2*params->n_particles);
  }
  else{
    fission_bank = init_bank(2*params->n_particles/n_threads);
  }
}
  master_fission_bank = init_bank(2*params->n_particles);
#else
  fission_bank = init_bank(2*params->n_particles);
#endif

  // Initialize source bank
  source_bank = init_bank(params->n_particles);

  // Sample source particles or load a source
  if(params->load_source == TRUE){
    load_source(source_bank);
    source_bank->n = params->n_particles;
  }
  else{
    for(i_p=0; i_p<params->n_particles; i_p++){
      sample_source_particle(&(source_bank->p[i_p]), g, &seed2);
      source_bank->n++;
    }
  }

  center_print("SIMULATION", 79);
  border_print();
  printf("%-15s %-15s %-15s %-15s\n", "BATCH", "ENTROPY", "KEFF", "MEAN KEFF");

  // Set initial seed before starting eigenvalue problem
  seed1_0 = seed1;
  seed2_0 = seed2;

  // Start time
  t1 = timer();

  // Loop over batches
  for(i_b=0; i_b<params->n_batches; i_b++){

    keff_batch = 0;

    // Write coordinates of particles in source bank
    if(params->write_bank == TRUE){
      write_bank(source_bank, fp, params->bank_file);
    }

    // Turn on tallying and increment index in active batches
    if(i_b >= params->n_batches - params->n_active){
      i_a++;
      if(params->tally == TRUE){
        t->tallies_on = TRUE;
      }
    }

    // Loop over generations
    for(i_g=0; i_g<params->n_generations; i_g++){

#pragma omp parallel default(none) private(i_p, p, seed1) shared(i_b, i_g, params, source_bank, g, m, t, seed1_0)
{
#pragma omp for
      // Loop over particles
      for(i_p=0; i_p<params->n_particles; i_p++){

	// Set seed for particle i_p by skipping ahead in the random number
	// sequence stride*(total particles simulated) numbers. This allows for
	// reproducibility of the particle history.
        seed1 = rn_skip(seed1_0, (i_b*params->n_generations + i_g)*params->n_particles + i_p);

        // Copy next particle into p
        copy_particle(&p, &(source_bank->p[i_p]));

        // Transport the next particle from source bank
	transport(&p, g, m, t, fission_bank, params, &seed1);
      }
}
      seed2 = rn_skip(seed2_0, i_b*params->n_generations + i_g);

// Merge fission banks from threads
#ifdef _OPENMP
      n_sites = 0;
#pragma omp parallel
{
#pragma omp for ordered
      for(i_t=0; i_t<n_threads; i_t++){
#pragma omp ordered
{
        memcpy(&(master_fission_bank->p[n_sites]), fission_bank->p, fission_bank->n*sizeof(Particle));
        n_sites += fission_bank->n;
}
      }
#pragma omp barrier
      // Copy shared fission bank sites into master thread's fission bank
      if(thread_id == 0){
        memcpy(fission_bank->p, master_fission_bank->p, n_sites*sizeof(Particle));
        fission_bank->n = n_sites;
      }
      else{
        fission_bank->n = 0;
      }
}
#endif

      // Calculate generation k_effective and accumulate batch k_effective
      keff_gen = (double) fission_bank->n / source_bank->n;
      keff_batch += keff_gen;

      // Sample new source particles from the particles that were added to the
      // fission bank during this generation
      synchronize_bank(source_bank, fission_bank, g, &seed2);

      // Calculate shannon entropy to assess source convergence
      H = shannon_entropy(g, source_bank);
      if(params->write_entropy == TRUE){
        write_entropy(H, fp, params->entropy_file);
      }

      // Write the source distribution
      if(params->write_source == TRUE){
        write_source(g, source_bank, params, fp, params->source_file);
      }
    }

    // Calculate k_effective
    keff_batch /= params->n_generations;
    if(i_a >= 0){
      keff[i_a] = keff_batch;
    }

    // Tallies for this realization
    if(t->tallies_on == TRUE){
      if(params->write_tally == TRUE){
        write_tally(t, fp, params->tally_file);
      }
      memset(t->flux, 0, t->n*t->n*t->n*sizeof(double));
    }

    // Calculate keff mean and standard deviation
    calculate_keff(keff, &keff_mean, &keff_std, i_a+1);

    // Status text
    if(i_a < 0){
      printf("%-15d %-15f %-15f\n", i_b+1, H, keff_batch);
    }
    else{
      printf("%-15d %-15f %-15f %f +/- %-15f\n", i_b+1, H, keff_batch, keff_mean, keff_std);
    }
  }

  // Stop time
  t2 = timer();

  printf("Simulation time: %f secs\n", t2-t1);

  // Write out keff
  if(params->write_keff == TRUE){
    write_keff(keff, params->n_active, fp, params->keff_file);
  }

  if(params->save_source == TRUE){
    save_source(source_bank);
  }

  // Free memory
#ifdef _OPENMP
#pragma omp parallel
{
  free_bank(fission_bank);
}
  free_bank(master_fission_bank);
#else
  free_bank(fission_bank);
#endif
  free_bank(source_bank);
  free_tally(t);
  free_material(m);
  free(g);
  free(keff);
  free(params);

  return 0;
}
