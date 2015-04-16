#include "simple_transport_header.h"

int main(int argc, char *argv[])
{
  //srand(time(NULL));

  int i_b;            // index over batches
  int i_g;            // index over generations
  unsigned long i_p;  // index over particles
  double t1, t2;      // timers
  double keff;        // k_effective
  double H;           // shannon entropy
  Parameters *params; // user defined parameters
  Geometry *g;
  Material *m;
  Bank *source_bank;
  Bank *fission_bank;

  // Get inputs
  params = set_default_params();
  parse_params("parameters", params);
  print_params(params);

  // Set up geometry
  g = init_geometry(params);

  // Set up material
  m = init_material(params);

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

    // Loop over generations
    for(i_g=0; i_g<params->n_generations; i_g++){

      // Loop over particles
      for(i_p=0; i_p<source_bank->n; i_p++){

        // Transport the next particle from source bank
        transport(&(source_bank->p[i_p]), g, m, fission_bank);
      }

      // Accumulate generation k_effective
      keff += (double) fission_bank->n / source_bank->n;

      // Sample new source particles from the particles that were added to the
      // fission bank during this generation
      synchronize_bank(source_bank, fission_bank, g);

      // Calculate shannon entropy to assess source convergence
      H = shannon_entropy(g, source_bank);
      printf("Shannon entropy: %f\n", H);
    }

    // Calculate k_effective
    keff /= params->n_generations;
    printf("keff: %f\n", keff);

  }

  // Stop time
  t2 = timer();

  printf("Simulation time: %f secs\n", t2-t1);

  // Free memory
  free_bank(source_bank);
  free_bank(fission_bank);
  free_material(m);
  free(g);

  return 0;
}
