#include "header.h"
#include "global.h"

int main(int argc, char *argv[])
{
  double t1, t2; // timers

  // Initialize the global variables used throughout the simulation
  init_problem(argc, argv);

  center_print("SIMULATION", 79);
  border_print();
  printf("%-15s %-15s %-15s %-15s\n", "BATCH", "ENTROPY", "KEFF", "MEAN KEFF");

  // Start time
  t1 = timer();

  run_eigenvalue();

  // Stop time
  t2 = timer();

  printf("Simulation time: %f secs\n", t2-t1);

  // Free memory
  free_problem();

  return 0;
}
