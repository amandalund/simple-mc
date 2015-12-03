#include "header.h"

int main(int argc, char *argv[])
{
  FILE *fp = NULL;    // file pointer for output
  double t1, t2;      // timers
  double *keff;       // effective multiplication factor
  Parameters *params; // user defined parameters
  Geometry *g;
  Material *m;
  Tally *t;
  Bank *source_bank;
  Bank *fission_bank;

  // Get inputs
  params = set_default_params();
  parse_params("parameters", params);
  read_CLI(argc, argv, params);
  print_params(params);

  // Set up output files
  init_output(params, fp);

  // Seed RNG
  srand(params->seed);

  // Set up array for keff
  keff = calloc(params->n_active, sizeof(double));

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

  center_print("SIMULATION", 79);
  border_print();
  printf("%-15s %-15s %-15s %-15s\n", "BATCH", "ENTROPY", "KEFF", "MEAN KEFF");

  // Start time
  t1 = timer();

  // Converge source (inactive batches)
  //converge_source(params, source_bank, fission_bank, g, m, t);
  ramp_up(params, source_bank, fission_bank, g, m, t);

  // Run eigenvalue problem (active batches)
  run_eigenvalue(params, source_bank, fission_bank, g, m, t, keff);

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
  free_bank(fission_bank);
  free_bank(source_bank);
  free_tally(t);
  free_material(m);
  free(g);
  free(keff);
  free(params);

  return 0;
}
