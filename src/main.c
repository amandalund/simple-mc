#include "simple_transport_header.h"

int main(int argc, char *argv[])
{
  int i_b;            // index over batches
  int i_g;            // index over generations
  unsigned long i_p;  // index over particles
  FILE *fp = NULL;    // file pointer for output
  double t1, t2;      // timers
  double keff;        // effective multiplication factor
  double H;           // shannon entropy
  Parameters *params; // user defined parameters
  Geometry *g;
  Material *m;
  Tally *t;
  Bank *source_bank;
  Bank *fission_bank;

  // Get inputs
  params = set_default_params();
  parse_params("parameters", params);
  print_params(params);

  // Setup output files
  init_output(params, fp);

  // Seed RNG
  srand(params->seed);

  // Set up geometry
  g = init_geometry(params);

  // Set up material
  m = init_material(params);

  // Set up tallies
  t = init_tally(params);

  // Initialize source bank
  source_bank = init_bank(params->n_particles);

  // Initialize fission bank
  fission_bank = init_bank(params->n_particles);

  // Sample source particles
  for(i_p=0; i_p<params->n_particles; i_p++){
    sample_source_particle(&(source_bank->p[i_p]), g);
    source_bank->n++;
  }

  // Start time
  t1 = timer();

  // Loop over batches
  for(i_b=0; i_b<params->n_batches; i_b++){

    keff = 0.0;

    // Turn on tallying
    if(params->tally == TRUE && i_b == params->n_batches - params->n_active){
      t->tallies_on = TRUE;
    }

    // Loop over generations
    for(i_g=0; i_g<params->n_generations; i_g++){

      // Loop over particles
      for(i_p=0; i_p<source_bank->n; i_p++){

        // Transport the next particle from source bank
        transport(&(source_bank->p[i_p]), g, m, t, fission_bank);
      }

      // Accumulate generation k_effective
      keff += (double) fission_bank->n / source_bank->n;

      // Sample new source particles from the particles that were added to the
      // fission bank during this generation
      synchronize_bank(source_bank, fission_bank, g);

      // Calculate shannon entropy to assess source convergence
      H = shannon_entropy(g, source_bank, params);
      if(params->write_entropy == TRUE){
        write_entropy(H, fp, params->entropy_file);
      }
    }

    // Calculate k_effective
    keff /= params->n_generations;
    if(params->write_keff == TRUE){
      write_keff(keff, fp, params->keff_file);
    }

    // Tallies for this realization
    if(t->tallies_on == TRUE){
      batch_tally(t, params);
      if(params->write_tally == TRUE){
        write_tally(t, fp, params->tally_file);
      }
    }

    // Status text
    printf("\rBatch %d/%d: keff = %f", i_b+1, params->n_batches, keff);
    fflush(stdout);
  }

  // Stop time
  t2 = timer();

  printf("\nSimulation time: %f secs\n", t2-t1);

  // Free memory
  free_bank(fission_bank);
  free_bank(source_bank);
  free_tally(t);
  free_material(m);
  free(g);
  free(params);

  return 0;
}
