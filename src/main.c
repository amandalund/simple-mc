#include "header.h"

int main(int argc, char *argv[])
{
  Parameters *parameters; // user defined parameters
  Geometry *geometry; // homogenous cube geometry
  Material *material; // problem material
  Bank *source_bank; // array for particle source sites
  Bank *fission_bank; // array for particle fission sites
  Tally *tally; // scalar flux tally
  double *keff; // effective multiplication factor
  double t1, t2; // timers

  // Get inputs: set parameters to default values, parse parameter file,
  // override with any command line inputs, and print parameters
  parameters = init_parameters();
  parse_parameters(parameters);
  read_CLI(argc, argv, parameters);
  print_parameters(parameters);

  // Set initial RNG seed
  set_initial_seed(parameters->seed);
  set_stream(STREAM_INIT);

  // Create files for writing results to
  init_output(parameters);

  // Set up geometry
  geometry = init_geometry(parameters);

  // Set up material
  material = init_material(parameters);

  // Set up tallies
  tally = init_tally(parameters);

  // Create source bank and initial source distribution
  source_bank = init_source_bank(parameters, geometry);

  // Create fission bank
  fission_bank = init_fission_bank(parameters);

  // Set up array for k effective
  keff = calloc(parameters->n_active, sizeof(double));

  center_print("SIMULATION", 79);
  border_print();
  printf("%-15s %-15s %-15s\n", "BATCH", "KEFF", "MEAN KEFF");

  // Start time
  t1 = timer();

  run_eigenvalue(parameters, geometry, material, source_bank, fission_bank, tally, keff);

  // Stop time
  t2 = timer();

  printf("Simulation time: %f secs\n", t2-t1);

  // Free memory
  free(keff);
  free_tally(tally);
  free_bank(fission_bank);
  free_bank(source_bank);
  free_material(material);
  free(geometry);
  free(parameters);

  return 0;
}
