#include "simple_transport_header.h"

// Reads command line inputs and applies options
void read_CLI(int argc, char *argv[], Input *input)
{
  int opt;
  while((opt = getopt(argc, argv, ":n:b:g:a:c:s:")) != -1){
    switch(opt){
      case 'n': input->n_particles = atol(optarg);   break;
      case 'b': input->n_batches = atoi(optarg);     break;
      case 'g': input->n_generations = atoi(optarg); break;
      case 'a': input->n_active = atoi(optarg);      break;
      case 'c': 
        if       (strcasecmp(optarg, "vacuum") == 0) input->bc = 0;
        else if (strcasecmp(optarg, "reflect") == 0) input->bc = 1;
        else if(strcasecmp(optarg, "periodic") == 0) input->bc = 2;
        else print_CLI_error();
        break;
      case 's': 
        if       (strcasecmp(optarg, "small") == 0) input->n_nuclides = 1;
        else if (strcasecmp(optarg, "medium") == 0) input->n_nuclides = 60;
        else if  (strcasecmp(optarg, "large") == 0) input->n_nuclides = 400;
        else print_CLI_error();
        break;
      default:  print_CLI_error();
    }
  } 

  // Validate inputs
  if((input->n_particles < 1) | (input->n_batches < 1 && input->n_generations < 1) |
     (input->n_batches < 0) | (input->n_generations < 0) | (input->n_active > input->n_batches)){
    print_CLI_error();
  }

  return;
}

// print error to screen, inform program options
void print_CLI_error(void)
{
  printf("Usage: ./simple_transport <options>\n");
  printf("Options include:\n");
  printf("  -n <particles>              Number of particles to simulate per fission source iteration\n");
  printf("  -b <batches>                Number of batches\n");
  printf("  -g <generations per batch>  Number of fission source iterations per batch\n");
  printf("  -a <active batches>         Number of batches which contribute to tallies\n");
  printf("  -c <boundary condition>     Boundary conditions (vacuum, reflect, periodic)\n");
  printf("  -s <material size>          Determines number of nuclides in the material (small, medium, large)\n");
  printf("See readme for full description of default run values\n");
  exit(1);
}

void print_inputs(Input *input)
{
  char *bc = NULL;
  if(input->bc == 0) bc = "Vacuum";
  else if(input->bc == 1) bc = "Reflective";
  else if(input->bc == 2) bc = "Periodic";
  border_print();
  center_print("INPUT SUMMARY", 79);
  border_print();
  printf("Number of particles:            "); fancy_int(input->n_particles);
  printf("Number of batches:              %d\n", input->n_batches);
  printf("Number of generations:          %d\n", input->n_generations);
  printf("Number of active batches:       %d\n", input->n_active);
  printf("Boundary conditions:            %s\n", bc);
  printf("Number of nuclides in material: %d\n", input->n_nuclides);
  border_print();
}

void border_print(void)
{
  printf( "=========================================="
     "======================================\n");
}

// Prints comma separated integers - for ease of reading
void fancy_int(long a)
{
  if(a < 1000)
    printf("%ld\n",a);
  else if(a >= 1000 && a < 1000000)
    printf("%ld,%03ld\n", a/1000, a % 1000);
  else if(a >= 1000000 && a < 1000000000)
    printf("%ld,%03ld,%03ld\n", a/1000000, (a % 1000000)/1000, a % 1000);
  else if(a >= 1000000000)
    printf("%ld,%03ld,%03ld,%03ld\n", a / 1000000000,
       (a % 1000000000)/1000000, (a % 1000000)/1000, a % 1000);
  else
    printf("%ld\n",a);
}

// Prints Section titles in center of 80 char terminal
void center_print(const char *s, int width)
{
  int length = strlen(s);
  int i;
  for (i=0; i<=(width-length)/2; i++){
    fputs(" ", stdout);
  }
  fputs(s, stdout);
  fputs("\n", stdout);
}

void print_particle(Particle *p)
{
  printf("alive: %-10dE: %-10.5fmu: %-10.5fphi: %-10.5fu: %-10.5f"
     "v: %-10.5fw: %-10.5fx: %-10.5fy: %-10.5fz: %-10.5f\n",
     p->alive, p->energy, p->mu, p->phi, p->u, p->v, p->w, p->x, p->y, p->z);

  return;
}

void print_bank(Bank *b)
{
  int i;
  Particle *p;

  for(i=0; i<b->n; i++){
    p = &(b->p[i]);
    print_particle(p);
  }

  return;
}
