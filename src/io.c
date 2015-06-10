#include "simple_transport_header.h"

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
    else
      printf("Unknown value '%s' in config file.\n", s);
  }

  fclose(fp);

  // Validate inputs
  if(params->write_tally == TRUE && params->tally_file == NULL)
    params->tally_file = "tally.dat";
  if(params->write_entropy == TRUE && params->entropy_file == NULL)
    params->entropy_file = "entropy.dat";
  if(params->write_keff == TRUE && params->keff_file == NULL)
    params->keff_file = "keff.dat";
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
  printf("Number of generations:          %d\n", params->n_generations);
  printf("Number of active batches:       %d\n", params->n_active);
  printf("Boundary conditions:            %s\n", bc);
  printf("Number of nuclides in material: %d\n", params->n_nuclides);
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

void write_keff(double keff, FILE *fp, char *filename)
{
  fp = fopen(filename, "a");
  fprintf(fp, "%.10f\n", keff);
  fclose(fp);

  return;
}
