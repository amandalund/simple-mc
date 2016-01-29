#include "header.h"

int main(int argc, char *argv[])
{
  FILE *fp = NULL;    // file pointer for output
  double t1, t2;      // timers
  double *keff;       // effective multiplication factor
  unsigned long i_p;  // index over particles
  Parameters *params; // user defined parameters
  Geometry *g;
  Material *m;
  Tally *t;
  Bank *source_bank;
  Bank *fission_bank;
  Queue *delay_queue;

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

  // Initialize delay queue
  delay_queue = init_queue(params->n_particles*params->lag);

  // Sample source particles or load a source
  if(params->load_source == TRUE){
    load_source(source_bank);
    source_bank->n = params->n_particles;
  }
  else{
    for(i_p=0; i_p<params->n_particles; i_p++){
      sample_source_particle(&(source_bank->p[i_p]), g, params);
      source_bank->n++;
    }
  }

  center_print("SIMULATION", 79);
  border_print();
  printf("%-15s %-15s %-15s %-15s\n", "BATCH", "ENTROPY", "KEFF", "MEAN KEFF");

  // Start time
  t1 = timer();

  // Converge source (inactive batches)
  printf("CONVERGING SOURCE...\n");
  converge_source(params, source_bank, fission_bank, g, m, t, delay_queue);

  // Build the delay bank before beginning active batches
  printf("BUILDING DELAY BANK...\n");
  build_delay_bank(params, source_bank, fission_bank, g, m, t, delay_queue);

  // Run eigenvalue problem (active batches)
  printf("BEGINNING ACTIVE CYCLES...\n");
  run_eigenvalue(params, source_bank, fission_bank, g, m, t, keff, delay_queue);

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
  free_queue(delay_queue);
  free_tally(t);
  free_material(m);
  free(g);
  free(keff);
  free(params);

  return 0;
}
