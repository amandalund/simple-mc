#include "header.h"
#include "global.h"

int main(int argc, char *argv[])
{
  int i_b; // index over batches
  int i_a = -1; // index over active batches
  int i_g; // index over generations
  unsigned long i_p; // index over particles
  FILE *fp = NULL; // file pointer for output
  double t1, t2; // timers
  double *keff; // effective multiplication factor
  double keff_gen = 1; // keff of generation
  double keff_batch; // keff of batch
  double keff_mean; // keff mean over active batches
  double keff_std; // keff standard deviation over active batches
  double H; // shannon entropy
  Particle p;

  // Initialize the global variables used throughout the simulation
  init_problem(argc, argv);

  // Set up array for keff
  keff = calloc(params->n_active, sizeof(double));

  center_print("SIMULATION", 79);
  border_print();
  printf("%-15s %-15s %-15s %-15s\n", "BATCH", "ENTROPY", "KEFF", "MEAN KEFF");

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

      // Set RNG stream for tracking
      set_stream(STREAM_TRACK);

#pragma omp parallel default(none) private(i_p, p) shared(i_b, i_g, params, source_bank, g, m, t)
{
#pragma omp for
      // Loop over particles
      for(i_p=0; i_p<params->n_particles; i_p++){

	// Set seed for particle i_p by skipping ahead in the random number
	// sequence stride*(total particles simulated) numbers. This allows for
	// reproducibility of the particle history.
        rn_skip((i_b*params->n_generations + i_g)*params->n_particles + i_p);

        // Copy next particle into p
        copy_particle(&p, &(source_bank->p[i_p]));

        // Transport the next particle from source bank
	transport(&p);
      }
}
      set_stream(STREAM_OTHER);
      rn_skip(i_b*params->n_generations + i_g);

#ifdef _OPENMP
      // Merge fission banks from threads
      merge_fission_banks();
#endif

      // Calculate generation k_effective and accumulate batch k_effective
      keff_gen = (double) fission_bank->n / source_bank->n;
      keff_batch += keff_gen;

      // Sample new source particles from the particles that were added to the
      // fission bank during this generation
      synchronize_bank();

      // Calculate shannon entropy to assess source convergence
      H = shannon_entropy(source_bank);
      if(params->write_entropy == TRUE){
        write_entropy(H, fp, params->entropy_file);
      }

      // Write the source distribution
      if(params->write_source == TRUE){
        write_source(source_bank, fp, params->source_file);
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
        write_tally(fp, params->tally_file);
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
  free_problem();
  free(keff);

  return 0;
}
