#include "header.h"

// Read in parameters from file
void parse_params(char *filename, Parameters *params)
{
  char line[256], *s;
  FILE *fp = fopen(filename, "r");

  while((s = fgets(line, sizeof(line), fp)) != NULL){

    if(line[0] == '#') continue;
    s = strtok(line, "=");
    if(s == NULL) continue;

    // Set parameters
    else if(strcmp(s, "particles") == 0){
      long long n_particles = atoll(strtok(NULL, "=\n"));
      if(n_particles < 1)
        print_error("Number of particles must be greater than 0");
      params->n_particles = n_particles;
    }
    else if(strcmp(s, "batches") == 0)
      params->n_batches = atoi(strtok(NULL, "=\n"));
    else if(strcmp(s, "generations") == 0)
      params->n_generations = atoi(strtok(NULL, "=\n"));
    else if(strcmp(s, "active") == 0)
      params->n_active = atoi(strtok(NULL, "=\n"));
    else if(strcmp(s, "nuclides") == 0)
      params->n_nuclides = atoi(strtok(NULL, "=\n"));
    else if(strcmp(s, "tally") == 0){
      s = strtok(NULL, "=\n");
      if(strcasecmp(s, "true") == 0)
        params->tally = TRUE;
      else if(strcasecmp(s, "false") == 0)
        params->tally = FALSE;
      else
        print_error("Invalid option for parameter 'tally': must be 'true' or 'false'");
    }
    else if(strcmp(s, "n_bins") == 0)
      params->n_bins = atoi(strtok(NULL, "=\n"));
    else if(strcmp(s, "seed") == 0)
      params->seed = atoi(strtok(NULL, "=\n"));
    else if(strcmp(s, "macro_xs_f") == 0)
      params->macro_xs_f = atof(strtok(NULL, "=\n"));
    else if(strcmp(s, "macro_xs_a") == 0)
      params->macro_xs_a = atof(strtok(NULL, "=\n"));
    else if(strcmp(s, "macro_xs_e") == 0)
      params->macro_xs_e = atof(strtok(NULL, "=\n"));
    else if(strcmp(s, "x") == 0)
      params->gx = atof(strtok(NULL, "=\n"));
    else if(strcmp(s, "y") == 0)
      params->gy = atof(strtok(NULL, "=\n"));
    else if(strcmp(s, "z") == 0)
      params->gz = atof(strtok(NULL, "=\n"));
    else if(strcmp(s, "bc") == 0){
      s = strtok(NULL, "=\n");
      if(strcasecmp(s, "vacuum") == 0)
        params->bc = 0;
      else if(strcasecmp(s, "reflective") == 0)
        params->bc = 1;
      else if(strcasecmp(s, "periodic") == 0)
        params->bc = 2;
      else
        print_error("Invalid boundary condition");
    }
    else if(strcmp(s, "write_tally") == 0){
      s = strtok(NULL, "=\n");
      if(strcasecmp(s, "true") == 0)
        params->write_tally = TRUE;
      else if(strcasecmp(s, "false") == 0)
        params->write_tally = FALSE;
      else
        print_error("Invalid option for parameter 'write_tally': must be 'true' or 'false'");
    }
    else if(strcmp(s, "write_entropy") == 0){
      s = strtok(NULL, "=\n");
      if(strcasecmp(s, "true") == 0)
        params->write_entropy = TRUE;
      else if(strcasecmp(s, "false") == 0)
        params->write_entropy = FALSE;
      else
        print_error("Invalid option for parameter 'write_entropy': must be 'true' or 'false'");
    }
    else if(strcmp(s, "write_keff") == 0){
      s = strtok(NULL, "=\n");
      if(strcasecmp(s, "true") == 0)
        params->write_keff = TRUE;
      else if(strcasecmp(s, "false") == 0)
        params->write_keff = FALSE;
      else
        print_error("Invalid option for parameter 'write_keff': must be 'true' or 'false'");
    }
    else if(strcmp(s, "write_bank") == 0){
      s = strtok(NULL, "=\n");
      if(strcasecmp(s, "true") == 0)
        params->write_bank = TRUE;
      else if(strcasecmp(s, "false") == 0)
        params->write_bank = FALSE;
      else
        print_error("Invalid option for parameter 'write_bank': must be 'true' or 'false'");
    }
    else if(strcmp(s, "tally_file") == 0){
      s = strtok(NULL, "=\n");
      params->tally_file = malloc(strlen(s)*sizeof(char)+1);
      strcpy(params->tally_file, s);
    }
    else if(strcmp(s, "entropy_file") == 0){
      s = strtok(NULL, "=\n");
      params->entropy_file = malloc(strlen(s)*sizeof(char)+1);
      strcpy(params->entropy_file, s);
    }
    else if(strcmp(s, "keff_file") == 0){
      s = strtok(NULL, "=\n");
      params->keff_file = malloc(strlen(s)*sizeof(char)+1);
      strcpy(params->keff_file, s);
    }
    else if(strcmp(s, "bank_file") == 0){
      s = strtok(NULL, "=\n");
      params->bank_file = malloc(strlen(s)*sizeof(char)+1);
      strcpy(params->bank_file, s);
    }
    else
      printf("Unknown value '%s' in config file.\n", s);
  }

  fclose(fp);

  return;
}

void read_CLI(int argc, char *argv[], Parameters *params)
{
  int i;
  char *arg;

  // Collect raw input
  for(i=1; i<argc; i++){
    arg = argv[i];

    // Number of particles (-p)
    if(strcmp(arg, "-p") == 0){
      if(++i < argc){
        long long n_particles = atoll(argv[i]);
        if(n_particles < 1)
          print_error("Number of particles must be greater than 0");
        params->n_particles = n_particles;
      }
      else print_error("Error reading command line input '-p'");
    }

    // Number of batches (-b)
    else if(strcmp(arg, "-b") == 0){
      if(++i < argc) params->n_batches = atoi(argv[i]);
      else print_error("Error reading command line input '-b'");
    }

    // Number of active batches (-a)
    else if(strcmp(arg, "-a") == 0){
      if(++i < argc) params->n_active = atoi(argv[i]);
      else print_error("Error reading command line input '-a'");
    }

    // Number of generations (-g)
    else if(strcmp(arg, "-g") == 0){
      if(++i < argc) params->n_generations = atoi(argv[i]);
      else print_error("Error reading command line input '-g'");
    }

    // Boundary conditions (-c)
    else if(strcmp(arg, "-c") == 0){
      if(++i < argc){
        if(strcasecmp(argv[i], "vacuum") == 0)
          params->bc = 0;
        else if(strcasecmp(argv[i], "reflective") == 0)
          params->bc = 1;
        else if(strcasecmp(argv[i], "periodic") == 0)
          params->bc = 2;
        else
          print_error("Invalid boundary condition");
      }
      else print_error("Error reading command line input '-c'");
    }

    // Number of nuclides in material (-n)
    else if(strcmp(arg, "-n") == 0){
      if(++i < argc) params->n_nuclides = atoi(argv[i]);
      else print_error("Error reading command line input '-n'");
    }

    // Whether to tally (-t)
    else if(strcmp(arg, "-t") == 0){
      if(++i < argc){
        if(strcasecmp(argv[i], "true") == 0)
          params->tally = TRUE;
        else if(strcasecmp(argv[i], "false") == 0)
          params->tally = FALSE;
        else
          print_error("Invalid option for parameter 'tally': must be 'true' or 'false'");
      }
      else print_error("Error reading command line input '-t'");
    }

    // Number of bins in each dimension (-m)
    else if(strcmp(arg, "-m") == 0){
      if(++i < argc) params->n_bins = atoi(argv[i]);
      else print_error("Error reading command line input '-m'");
    }

    // RNG seed (-s)
    else if(strcmp(arg, "-s") == 0){
      if(++i < argc) params->seed = atoi(argv[i]);
      else print_error("Error reading command line input '-s'");
    }

    // Number of batches (-b)
    else if(strcmp(arg, "-b") == 0){
      if(++i < argc) params->n_batches = atoi(argv[i]);
      else print_error("Error reading command line input '-b'");
    }

    // Absorption macro xs (-d)
    else if(strcmp(arg, "-d") == 0){
      if(++i < argc) params->macro_xs_a = atof(argv[i]);
      else print_error("Error reading command line input '-d'");
    }

    // Elastic macro xs (-e)
    else if(strcmp(arg, "-e") == 0){
      if(++i < argc) params->macro_xs_e = atof(argv[i]);
      else print_error("Error reading command line input '-e'");
    }

    // Fission macro xs (-f)
    else if(strcmp(arg, "-f") == 0){
      if(++i < argc) params->macro_xs_f = atof(argv[i]);
      else print_error("Error reading command line input '-f'");
    }

    // Geometry size in x (-x)
    else if(strcmp(arg, "-x") == 0){
      if(++i < argc) params->gx = atof(argv[i]);
      else print_error("Error reading command line input '-x'");
    }

    // Geometry size in y (-y)
    else if(strcmp(arg, "-y") == 0){
      if(++i < argc) params->gy = atof(argv[i]);
      else print_error("Error reading command line input '-y'");
    }

    // Geometry size in z (-z)
    else if(strcmp(arg, "-z") == 0){
      if(++i < argc) params->gz = atof(argv[i]);
      else print_error("Error reading command line input '-z'");
    }

    // Whether to output tally (-i)
    else if(strcmp(arg, "-i") == 0){
      if(++i < argc){
        if(strcasecmp(argv[i], "true") == 0)
          params->write_tally = TRUE;
        else if(strcasecmp(argv[i], "false") == 0)
          params->write_tally = FALSE;
        else
          print_error("Invalid option for parameter 'write_tally': must be 'true' or 'false'");
      }
      else print_error("Error reading command line input '-i'");
    }

    // Whether to output shannon entropy (-j)
    else if(strcmp(arg, "-j") == 0){
      if(++i < argc){
        if(strcasecmp(argv[i], "true") == 0)
          params->write_entropy = TRUE;
        else if(strcasecmp(argv[i], "false") == 0)
          params->write_entropy = FALSE;
        else
          print_error("Invalid option for parameter 'write_entropy': must be 'true' or 'false'");
      }
      else print_error("Error reading command line input '-j'");
    }

    // Whether to output keff (-k)
    else if(strcmp(arg, "-k") == 0){
      if(++i < argc){
        if(strcasecmp(argv[i], "true") == 0)
          params->write_keff = TRUE;
        else if(strcasecmp(argv[i], "false") == 0)
          params->write_keff = FALSE;
        else
          print_error("Invalid option for parameter 'write_keff': must be 'true' or 'false'");
      }
      else print_error("Error reading command line input '-k'");
    }

    // Whether to output particle bank (-q)
    else if(strcmp(arg, "-q") == 0){
      if(++i < argc){
        if(strcasecmp(argv[i], "true") == 0)
          params->write_bank = TRUE;
        else if(strcasecmp(argv[i], "false") == 0)
          params->write_bank = FALSE;
        else
          print_error("Invalid option for parameter 'write_bank': must be 'true' or 'false'");
      }
      else print_error("Error reading command line input '-q'");
    }

    // Path to write tallies to (-u)
    else if(strcmp(arg, "-u") == 0){
      if(++i < argc){
        if(params->tally_file != NULL) free(params->tally_file);
        params->tally_file = malloc(strlen(argv[i])*sizeof(char)+1);
        strcpy(params->tally_file, argv[i]);
      }
      else print_error("Error reading command line input '-u'");
    }

    // Path to write shannon entropy to (-v)
    else if(strcmp(arg, "-v") == 0){
      if(++i < argc){
        if(params->entropy_file != NULL) free(params->entropy_file);
        params->entropy_file = malloc(strlen(argv[i])*sizeof(char)+1);
        strcpy(params->entropy_file, argv[i]);
      }
      else print_error("Error reading command line input '-v'");
    }

    // Path to write keff to (-w)
    else if(strcmp(arg, "-w") == 0){
      if(++i < argc){
        if(params->keff_file != NULL) free(params->keff_file);
        params->keff_file = malloc(strlen(argv[i])*sizeof(char)+1);
        strcpy(params->keff_file, argv[i]);
      }
      else print_error("Error reading command line input '-w'");
    }

    // Path to write bank to (-r)
    else if(strcmp(arg, "-r") == 0){
      if(++i < argc){
        if(params->bank_file != NULL) free(params->bank_file);
        params->bank_file = malloc(strlen(argv[i])*sizeof(char)+1);
        strcpy(params->bank_file, argv[i]);
      }
      else print_error("Error reading command line input '-r'");
    }

    else print_error("Error reading command line input");
  }

  // Validate Inputs
  if(params->write_tally == TRUE && params->tally_file == NULL)
    params->tally_file = "tally.dat";
  if(params->write_entropy == TRUE && params->entropy_file == NULL)
    params->entropy_file = "entropy.dat";
  if(params->write_keff == TRUE && params->keff_file == NULL)
    params->keff_file = "keff.dat";
  if(params->write_bank == TRUE && params->bank_file == NULL)
    params->bank_file = "bank.dat";
  if(params->n_batches < 1 && params->n_generations < 1)
    print_error("Must have at least one batch or one generation");
  if(params->n_batches < 0)
    print_error("Number of batches cannot be negative");
  if(params->n_generations < 0)
    print_error("Number of generations cannot be negative");
  if(params->n_active > params->n_batches)
    print_error("Number of active batches cannot be greater than number of batches");
  if(params->n_bins < 0)
    print_error("Number of bins cannot be negative");
  if(params->gx <= 0 || params->gy <= 0 || params->gz <= 0)
    print_error("Length of domain must be positive in x, y, and z dimension");
  if(params->macro_xs_f < 0 || params->macro_xs_a < 0 || params->macro_xs_e < 0)
    print_error("Macroscopic cross section values cannot be negative");

  return;
}

void print_error(char *message)
{
  printf("ERROR: %s\n", message);
  exit(1);
}

void print_params(Parameters *params)
{
  char *bc = NULL;
  if(params->bc == 0) bc = "Vacuum";
  else if(params->bc == 1) bc = "Reflective";
  else if(params->bc == 2) bc = "Periodic";
  border_print();
  center_print("INPUT SUMMARY", 79);
  border_print();
  printf("Number of particles:            "); fancy_int(params->n_particles);
  printf("Number of batches:              %d\n", params->n_batches);
  printf("Number of active batches:       %d\n", params->n_active);
  printf("Number of generations:          %d\n", params->n_generations);
  printf("Boundary conditions:            %s\n", bc);
  printf("Number of nuclides in material: %d\n", params->n_nuclides);
  printf("RNG seed:                       %d\n", params->seed);
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

void write_tally(Tally *t, FILE *fp, char *filename)
{
  int i, j, k;

  fp = fopen(filename, "a");

  for(i=0; i<t->n; i++){
    for(j=0; j<t->n; j++){
      for(k=0; k<t->n; k++){
        fprintf(fp, "%e ", t->mean[i + t->n*j + t->n*t->n*k]);
      }
      fprintf(fp, "\n");
    }
  }

  fclose(fp);

  return;
}

void write_entropy(double H, FILE *fp, char *filename)
{
  fp = fopen(filename, "a");
  fprintf(fp, "%.10f\n", H);
  fclose(fp);

  return;
}

void write_keff(double *keff, int n, FILE *fp, char *filename)
{
  int i;
  fp = fopen(filename, "a");

  for(i=0; i<n; i++){
    fprintf(fp, "%.10f\n", keff[i]);
  }

  fclose(fp);

  return;
}

void write_bank(Bank *b, FILE *fp, char *filename)
{
  int i;

  fp = fopen(filename, "a");

  for(i=0; i<b->n; i++){
    fprintf(fp, "%.10f %.10f %.10f ", b->p[i].x, b->p[i].y, b->p[i].z);
  }
  fprintf(fp, "\n");

  fclose(fp);
  return;
}
