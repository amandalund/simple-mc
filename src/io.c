#include "header.h"

// Read in parameters from file
void parse_parameters(Parameters *parameters)
{
  char line[256], *s;
  FILE *fp = fopen("parameters", "r");

  while((s = fgets(line, sizeof(line), fp)) != NULL){

    if(line[0] == '#') continue;
    if(line[0] == '\n') continue;
    s = strtok(line, "=");
    if(s == NULL) continue;

    // Number of particles
    else if(strcmp(s, "particles") == 0){
      long long n_particles = atoll(strtok(NULL, "=\n"));
      if(n_particles < 1)
        print_error("Number of particles must be greater than 0");
      parameters->n_particles = n_particles;
    }

    // Number of batches
    else if(strcmp(s, "batches") == 0){
      parameters->n_batches = atoi(strtok(NULL, "=\n"));
    }

    // Number of generations
    else if(strcmp(s, "generations") == 0){
      parameters->n_generations = atoi(strtok(NULL, "=\n"));
    }

    // Number of active batches
    else if(strcmp(s, "active") == 0){
      parameters->n_active = atoi(strtok(NULL, "=\n"));
    }

    // Number of nuclides in material
    else if(strcmp(s, "nuclides") == 0){
      parameters->n_nuclides = atoi(strtok(NULL, "=\n"));
    }

    // Whether to tally
    else if(strcmp(s, "tally") == 0){
      s = strtok(NULL, "=\n");
      if(strcasecmp(s, "true") == 0)
        parameters->tally = TRUE;
      else if(strcasecmp(s, "false") == 0)
        parameters->tally = FALSE;
      else
        print_error("Invalid option for parameter 'tally': must be 'true' or 'false'");
    }

    // Number of bins in each dimension
    else if(strcmp(s, "bins") == 0){
      parameters->n_bins = atoi(strtok(NULL, "=\n"));
    }

    // RNG seed
    else if(strcmp(s, "seed") == 0){
      parameters->seed = atol(strtok(NULL, "=\n"));
    }

    // Average number of fission neutrons produced
    else if(strcmp(s, "nu") == 0){
      parameters->nu = atof(strtok(NULL, "=\n"));
    }

    // Fission macro xs
    else if(strcmp(s, "xs_f") == 0){
      parameters->xs_f = atof(strtok(NULL, "=\n"));
    }

    // Absorption macro xs
    else if(strcmp(s, "xs_a") == 0){
      parameters->xs_a = atof(strtok(NULL, "=\n"));
    }

    // Scattering macro xs
    else if(strcmp(s, "xs_s") == 0){
      parameters->xs_s = atof(strtok(NULL, "=\n"));
    }

    // Geometry size in x
    else if(strcmp(s, "x") == 0){
      parameters->gx = atof(strtok(NULL, "=\n"));
    }

    // Geometry size in y
    else if(strcmp(s, "y") == 0){
      parameters->gy = atof(strtok(NULL, "=\n"));
    }

    // Geometry size in z
    else if(strcmp(s, "z") == 0){
      parameters->gz = atof(strtok(NULL, "=\n"));
    }

    // Boundary conditions
    else if(strcmp(s, "bc") == 0){
      s = strtok(NULL, "=\n");
      if(strcasecmp(s, "vacuum") == 0)
        parameters->bc = 0;
      else if(strcasecmp(s, "reflective") == 0)
        parameters->bc = 1;
      else if(strcasecmp(s, "periodic") == 0)
        parameters->bc = 2;
      else
        print_error("Invalid boundary condition");
    }

    // Whether to output tally
    else if(strcmp(s, "write_tally") == 0){
      s = strtok(NULL, "=\n");
      if(strcasecmp(s, "true") == 0)
        parameters->write_tally = TRUE;
      else if(strcasecmp(s, "false") == 0)
        parameters->write_tally = FALSE;
      else
        print_error("Invalid option for parameter 'write_tally': must be 'true' or 'false'");
    }

    // Whether to output keff
    else if(strcmp(s, "write_keff") == 0){
      s = strtok(NULL, "=\n");
      if(strcasecmp(s, "true") == 0)
        parameters->write_keff = TRUE;
      else if(strcasecmp(s, "false") == 0)
        parameters->write_keff = FALSE;
      else
        print_error("Invalid option for parameter 'write_keff': must be 'true' or 'false'");
    }

    // Path to write tallies to
    else if(strcmp(s, "tally_file") == 0){
      s = strtok(NULL, "=\n");
      parameters->tally_file = malloc(strlen(s)*sizeof(char)+1);
      strcpy(parameters->tally_file, s);
    }

    // Path to write keff to
    else if(strcmp(s, "keff_file") == 0){
      s = strtok(NULL, "=\n");
      parameters->keff_file = malloc(strlen(s)*sizeof(char)+1);
      strcpy(parameters->keff_file, s);
    }

    // Unknown config file option
    else print_error("Unknown option in config file.");
  }

  fclose(fp);

  return;
}

void read_CLI(int argc, char *argv[], Parameters *parameters)
{
  int i;
  char *arg;

  // Collect raw input
  for(i=1; i<argc; i++){
    arg = argv[i];

    // Number of particles (-particles)
    if(strcmp(arg, "-particles") == 0){
      if(++i < argc){
        long long n_particles = atoll(argv[i]);
        if(n_particles < 1)
          print_error("Number of particles must be greater than 0");
        parameters->n_particles = n_particles;
      }
      else print_error("Error reading command line input '-particles'");
    }

    // Number of batches (-batches)
    else if(strcmp(arg, "-batches") == 0){
      if(++i < argc) parameters->n_batches = atoi(argv[i]);
      else print_error("Error reading command line input '-batches'");
    }

    // Number of active batches (-active)
    else if(strcmp(arg, "-active") == 0){
      if(++i < argc) parameters->n_active = atoi(argv[i]);
      else print_error("Error reading command line input '-active'");
    }

    // Number of generations (-generations)
    else if(strcmp(arg, "-generations") == 0){
      if(++i < argc) parameters->n_generations = atoi(argv[i]);
      else print_error("Error reading command line input '-generations'");
    }

    // Boundary conditions (-bc)
    else if(strcmp(arg, "-bc") == 0){
      if(++i < argc){
        if(strcasecmp(argv[i], "vacuum") == 0)
          parameters->bc = 0;
        else if(strcasecmp(argv[i], "reflective") == 0)
          parameters->bc = 1;
        else if(strcasecmp(argv[i], "periodic") == 0)
          parameters->bc = 2;
        else
          print_error("Invalid boundary condition");
      }
      else print_error("Error reading command line input '-bc'");
    }

    // Number of nuclides in material (-nuclides)
    else if(strcmp(arg, "-nuclides") == 0){
      if(++i < argc) parameters->n_nuclides = atoi(argv[i]);
      else print_error("Error reading command line input '-nuclides'");
    }

    // Whether to tally (-tally)
    else if(strcmp(arg, "-tally") == 0){
      if(++i < argc){
        if(strcasecmp(argv[i], "true") == 0)
          parameters->tally = TRUE;
        else if(strcasecmp(argv[i], "false") == 0)
          parameters->tally = FALSE;
        else
          print_error("Invalid option for parameter 'tally': must be 'true' or 'false'");
      }
      else print_error("Error reading command line input '-tally'");
    }

    // Number of bins in each dimension (-bins)
    else if(strcmp(arg, "-bins") == 0){
      if(++i < argc) parameters->n_bins = atoi(argv[i]);
      else print_error("Error reading command line input '-bins'");
    }

    // RNG seed (-seed)
    else if(strcmp(arg, "-seed") == 0){
      if(++i < argc) parameters->seed = atol(argv[i]);
      else print_error("Error reading command line input '-seed'");
    }

    // Average number of fission neutrons produced (-nu)
    else if(strcmp(arg, "-nu") == 0){
      if(++i < argc) parameters->nu = atof(argv[i]);
      else print_error("Error reading command line input '-nu'");
    }

    // Absorption macro xs (-xs_a)
    else if(strcmp(arg, "-xs_a") == 0){
      if(++i < argc) parameters->xs_a = atof(argv[i]);
      else print_error("Error reading command line input '-xs_a'");
    }

    // Scattering macro xs (-xs_s)
    else if(strcmp(arg, "-xs_s") == 0){
      if(++i < argc) parameters->xs_s = atof(argv[i]);
      else print_error("Error reading command line input '-xs_s'");
    }

    // Fission macro xs (-xs_f)
    else if(strcmp(arg, "-xs_f") == 0){
      if(++i < argc) parameters->xs_f = atof(argv[i]);
      else print_error("Error reading command line input '-xs_f'");
    }

    // Geometry size in x (-x)
    else if(strcmp(arg, "-x") == 0){
      if(++i < argc) parameters->gx = atof(argv[i]);
      else print_error("Error reading command line input '-x'");
    }

    // Geometry size in y (-y)
    else if(strcmp(arg, "-y") == 0){
      if(++i < argc) parameters->gy = atof(argv[i]);
      else print_error("Error reading command line input '-y'");
    }

    // Geometry size in z (-z)
    else if(strcmp(arg, "-z") == 0){
      if(++i < argc) parameters->gz = atof(argv[i]);
      else print_error("Error reading command line input '-z'");
    }

    // Whether to output tally (-write_tally)
    else if(strcmp(arg, "-write_tally") == 0){
      if(++i < argc){
        if(strcasecmp(argv[i], "true") == 0)
          parameters->write_tally = TRUE;
        else if(strcasecmp(argv[i], "false") == 0)
          parameters->write_tally = FALSE;
        else
          print_error("Invalid option for parameter 'write_tally': must be 'true' or 'false'");
      }
      else print_error("Error reading command line input '-write_tally'");
    }

    // Whether to output keff (-write_keff)
    else if(strcmp(arg, "-write_keff") == 0){
      if(++i < argc){
        if(strcasecmp(argv[i], "true") == 0)
          parameters->write_keff = TRUE;
        else if(strcasecmp(argv[i], "false") == 0)
          parameters->write_keff = FALSE;
        else
          print_error("Invalid option for parameter 'write_keff': must be 'true' or 'false'");
      }
      else print_error("Error reading command line input '-write_keff'");
    }

    // Path to write tallies to (-tally_file)
    else if(strcmp(arg, "-tally_file") == 0){
      if(++i < argc){
        if(parameters->tally_file != NULL) free(parameters->tally_file);
        parameters->tally_file = malloc(strlen(argv[i])*sizeof(char)+1);
        strcpy(parameters->tally_file, argv[i]);
      }
      else print_error("Error reading command line input '-tally_file'");
    }

    // Path to write keff to (-keff_file)
    else if(strcmp(arg, "-keff_file") == 0){
      if(++i < argc){
        if(parameters->keff_file != NULL) free(parameters->keff_file);
        parameters->keff_file = malloc(strlen(argv[i])*sizeof(char)+1);
        strcpy(parameters->keff_file, argv[i]);
      }
      else print_error("Error reading command line input '-keff_file'");
    }

    // Unknown command line option
    else print_error("Error reading command line input");
  }

  // Validate Inputs
  if(parameters->write_tally == TRUE && parameters->tally_file == NULL)
    parameters->tally_file = "tally.dat";
  if(parameters->write_keff == TRUE && parameters->keff_file == NULL)
    parameters->keff_file = "keff.dat";
  if(parameters->n_batches < 1 && parameters->n_generations < 1)
    print_error("Must have at least one batch or one generation");
  if(parameters->n_batches < 0)
    print_error("Number of batches cannot be negative");
  if(parameters->n_generations < 0)
    print_error("Number of generations cannot be negative");
  if(parameters->n_active > parameters->n_batches)
    print_error("Number of active batches cannot be greater than number of batches");
  if(parameters->n_bins < 0)
    print_error("Number of bins cannot be negative");
  if(parameters->nu < 0)
    print_error("Average number of fission neutrons produced cannot be negative");
  if(parameters->gx <= 0 || parameters->gy <= 0 || parameters->gz <= 0)
    print_error("Length of domain must be positive in x, y, and z dimension");
  if(parameters->xs_f < 0 || parameters->xs_a < 0 || parameters->xs_s < 0)
    print_error("Macroscopic cross section values cannot be negative");

  return;
}

void print_parameters(Parameters *parameters)
{
  char *bc = NULL;
  if(parameters->bc == 0) bc = "Vacuum";
  else if(parameters->bc == 1) bc = "Reflective";
  else if(parameters->bc == 2) bc = "Periodic";
  border_print();
  center_print("INPUT SUMMARY", 79);
  border_print();
  printf("Number of particles:            "); fancy_int(parameters->n_particles);
  printf("Number of batches:              %d\n", parameters->n_batches);
  printf("Number of active batches:       %d\n", parameters->n_active);
  printf("Number of generations:          %d\n", parameters->n_generations);
  printf("Boundary conditions:            %s\n", bc);
  printf("Number of nuclides in material: %d\n", parameters->n_nuclides);
  printf("RNG seed:                       %llu\n", parameters->seed);
  border_print();
}

void print_error(char *message)
{
  printf("ERROR: %s\n", message);
  exit(1);
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

void print_status(int i_a, int i_b, double keff_batch, double keff_mean, double keff_std)
{
  if(i_a < 0){
    printf("%-15d %-15f\n", i_b+1, keff_batch);
  }
  else{
    printf("%-15d %-15f %f +/- %-15f\n", i_b+1, keff_batch, keff_mean, keff_std);
  }

  return;
}

void init_output(Parameters *parameters)
{
  FILE *fp = NULL; // file pointer for output

  // Set up file to output tallies
  if(parameters->write_tally == TRUE){
    fp = fopen(parameters->tally_file, "w");
    fclose(fp);
  }

  // Set up file to output keff
  if(parameters->write_keff == TRUE){
    fp = fopen(parameters->keff_file, "w");
    fclose(fp);
  }

  return;
}

void write_tally(Tally *t, char *filename)
{
  int i, j, k;
  FILE *fp;

  fp = fopen(filename, "a");

  for(i=0; i<t->n; i++){
    for(j=0; j<t->n; j++){
      for(k=0; k<t->n; k++){
        fprintf(fp, "%e ", t->flux[i + t->n*j + t->n*t->n*k]);
      }
      fprintf(fp, "\n");
    }
  }

  fclose(fp);

  return;
}

void write_keff(double *keff, int n, char *filename)
{
  int i;
  FILE *fp;

  fp = fopen(filename, "a");

  for(i=0; i<n; i++){
    fprintf(fp, "%.10f\n", keff[i]);
  }

  fclose(fp);

  return;
}

