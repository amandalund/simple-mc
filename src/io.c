#include "simple_mc.h"

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

    // Number of bins in each dimension for shannon entropy
    else if(strcmp(s, "entropy_bins") == 0){
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

    // Domain length in x
    else if(strcmp(s, "Lx") == 0){
      parameters->Lx = atof(strtok(NULL, "=\n"));
    }

    // Domain length in y
    else if(strcmp(s, "Ly") == 0){
      parameters->Ly = atof(strtok(NULL, "=\n"));
    }

    // Domain length in z
    else if(strcmp(s, "Lz") == 0){
      parameters->Lz = atof(strtok(NULL, "=\n"));
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

    // Whether to load source
    else if(strcmp(s, "load_source") == 0){
      s = strtok(NULL, "=\n");
      if(strcasecmp(s, "true") == 0)
        parameters->load_source = TRUE;
      else if(strcasecmp(s, "false") == 0)
        parameters->load_source = FALSE;
      else
        print_error("Invalid option for parameter 'load_source': must be 'true' or 'false'");
    }

    // Whether to save source
    else if(strcmp(s, "save_source") == 0){
      s = strtok(NULL, "=\n");
      if(strcasecmp(s, "true") == 0)
        parameters->save_source = TRUE;
      else if(strcasecmp(s, "false") == 0)
        parameters->save_source = FALSE;
      else
        print_error("Invalid option for parameter 'save_source': must be 'true' or 'false'");
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

    // Whether to output shannon entropy
    else if(strcmp(s, "write_entropy") == 0){
      s = strtok(NULL, "=\n");
      if(strcasecmp(s, "true") == 0)
        parameters->write_entropy = TRUE;
      else if(strcasecmp(s, "false") == 0)
        parameters->write_entropy = FALSE;
      else
        print_error("Invalid option for parameter 'write_entropy': must be 'true' or 'false'");
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

    // Whether to output particle bank
    else if(strcmp(s, "write_bank") == 0){
      s = strtok(NULL, "=\n");
      if(strcasecmp(s, "true") == 0)
        parameters->write_bank = TRUE;
      else if(strcasecmp(s, "false") == 0)
        parameters->write_bank = FALSE;
      else
        print_error("Invalid option for parameter 'write_bank': must be 'true' or 'false'");
    }

    // Whether to output source distribution
    else if(strcmp(s, "write_source") == 0){
      s = strtok(NULL, "=\n");
      if(strcasecmp(s, "true") == 0)
        parameters->write_source = TRUE;
      else if(strcasecmp(s, "false") == 0)
        parameters->write_source = FALSE;
      else
        print_error("Invalid option for parameter 'write_source': must be 'true' or 'false'");
    }

    // Whether to output histories
    else if(strcmp(s, "write_histories") == 0){
      s = strtok(NULL, "=\n");
      if(strcasecmp(s, "true") == 0)
        parameters->write_histories = TRUE;
      else if(strcasecmp(s, "false") == 0)
        parameters->write_histories = FALSE;
      else
        print_error("Invalid option for parameter 'write_histories': must be 'true' or 'false'");
    }

    // Path to write tallies to
    else if(strcmp(s, "tally_file") == 0){
      s = strtok(NULL, "=\n");
      parameters->tally_file = malloc(strlen(s)*sizeof(char)+1);
      strcpy(parameters->tally_file, s);
    }

    // Path to write shannon entropy to
    else if(strcmp(s, "entropy_file") == 0){
      s = strtok(NULL, "=\n");
      parameters->entropy_file = malloc(strlen(s)*sizeof(char)+1);
      strcpy(parameters->entropy_file, s);
    }

    // Path to write keff to
    else if(strcmp(s, "keff_file") == 0){
      s = strtok(NULL, "=\n");
      parameters->keff_file = malloc(strlen(s)*sizeof(char)+1);
      strcpy(parameters->keff_file, s);
    }

    // Path to write bank to
    else if(strcmp(s, "bank_file") == 0){
      s = strtok(NULL, "=\n");
      parameters->bank_file = malloc(strlen(s)*sizeof(char)+1);
      strcpy(parameters->bank_file, s);
    }

    // Path to write source distribution to
    else if(strcmp(s, "source_file") == 0){
      s = strtok(NULL, "=\n");
      parameters->source_file = malloc(strlen(s)*sizeof(char)+1);
      strcpy(parameters->source_file, s);
    }

    // Path to write histories to
    else if(strcmp(s, "histories_file") == 0){
      s = strtok(NULL, "=\n");
      parameters->histories_file = malloc(strlen(s)*sizeof(char)+1);
      strcpy(parameters->histories_file, s);
    }

    // Convergence method
    else if(strcmp(s, "ramp_up") == 0){
      s = strtok(NULL, "=\n");
      if(strcasecmp(s, "false") == 0)
        parameters->ramp_up = FALSE;
      else if(strcasecmp(s, "true") == 0)
        parameters->ramp_up = TRUE;
      else
        print_error("Invalid convergence method: options for 'ramp_up' are <true>, <false>");
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

    // Number of bins in each dimension for shannon entropy  (-entropy_bins)
    else if(strcmp(arg, "-entropy_bins") == 0){
      if(++i < argc) parameters->entropy_bins = atoi(argv[i]);
      else print_error("Error reading command line input '-entropy_bins'");
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

    // Domain length in x (-Lx)
    else if(strcmp(arg, "-Lx") == 0){
      if(++i < argc) parameters->Lx = atof(argv[i]);
      else print_error("Error reading command line input '-Lx'");
    }

    // Domain length in y (-Ly)
    else if(strcmp(arg, "-Ly") == 0){
      if(++i < argc) parameters->Ly = atof(argv[i]);
      else print_error("Error reading command line input '-Ly'");
    }

    // Geometry size in z (-Lz)
    else if(strcmp(arg, "-Lz") == 0){
      if(++i < argc) parameters->Lz = atof(argv[i]);
      else print_error("Error reading command line input '-Lz'");
    }

    // Whether to load source (-load_source)
    else if(strcmp(arg, "-load_source") == 0){
      if(++i < argc){
        if(strcasecmp(argv[i], "true") == 0)
          parameters->load_source = TRUE;
        else if(strcasecmp(argv[i], "false") == 0)
          parameters->load_source = FALSE;
        else
          print_error("Invalid option for parameter 'load_source': must be 'true' or 'false'");
      }
      else print_error("Error reading command line input '-load_source'");
    }

    // Whether to save source (-save_source)
    else if(strcmp(arg, "-save_source") == 0){
      if(++i < argc){
        if(strcasecmp(argv[i], "true") == 0)
          parameters->save_source = TRUE;
        else if(strcasecmp(argv[i], "false") == 0)
          parameters->save_source = FALSE;
        else
          print_error("Invalid option for parameter 'save_source': must be 'true' or 'false'");
      }
      else print_error("Error reading command line input '-save_source'");
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

    // Whether to output shannon entropy (-write_entropy)
    else if(strcmp(arg, "-write_entropy") == 0){
      if(++i < argc){
        if(strcasecmp(argv[i], "true") == 0)
          parameters->write_entropy = TRUE;
        else if(strcasecmp(argv[i], "false") == 0)
          parameters->write_entropy = FALSE;
        else
          print_error("Invalid option for parameter 'write_entropy': must be 'true' or 'false'");
      }
      else print_error("Error reading command line input '-write_entropy'");
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

    // Whether to output particle bank (-write_bank)
    else if(strcmp(arg, "-write_bank") == 0){
      if(++i < argc){
        if(strcasecmp(argv[i], "true") == 0)
          parameters->write_bank = TRUE;
        else if(strcasecmp(argv[i], "false") == 0)
          parameters->write_bank = FALSE;
        else
          print_error("Invalid option for parameter 'write_bank': must be 'true' or 'false'");
      }
      else print_error("Error reading command line input '-write_bank'");
    }

    // Whether to output source distribution (-write_source)
    else if(strcmp(arg, "-write_source") == 0){
      if(++i < argc){
        if(strcasecmp(argv[i], "true") == 0)
          parameters->write_source = TRUE;
        else if(strcasecmp(argv[i], "false") == 0)
          parameters->write_source = FALSE;
        else
          print_error("Invalid option for parameter 'write_source': must be 'true' or 'false'");
      }
      else print_error("Error reading command line input '-write_source'");
    }

    // Whether to output histories (-write_histories)
    else if(strcmp(arg, "-write_histories") == 0){
      if(++i < argc){
        if(strcasecmp(argv[i], "true") == 0)
          parameters->write_histories = TRUE;
        else if(strcasecmp(argv[i], "false") == 0)
          parameters->write_histories = FALSE;
        else
          print_error("Invalid option for parameter 'write_histories': must be 'true' or 'false'");
      }
      else print_error("Error reading command line input '-write_source'");
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

    // Path to write shannon entropy to (-entropy_file)
    else if(strcmp(arg, "-entropy_file") == 0){
      if(++i < argc){
        if(parameters->entropy_file != NULL) free(parameters->entropy_file);
        parameters->entropy_file = malloc(strlen(argv[i])*sizeof(char)+1);
        strcpy(parameters->entropy_file, argv[i]);
      }
      else print_error("Error reading command line input '-entropy_file'");
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

    // Path to write bank to (-bank_file)
    else if(strcmp(arg, "-bank_file") == 0){
      if(++i < argc){
        if(parameters->bank_file != NULL) free(parameters->bank_file);
        parameters->bank_file = malloc(strlen(argv[i])*sizeof(char)+1);
        strcpy(parameters->bank_file, argv[i]);
      }
      else print_error("Error reading command line input '-bank_file'");
    }

    // Path to write source distribution to (-source_file)
    else if(strcmp(arg, "-source_file") == 0){
      if(++i < argc){
        if(parameters->source_file != NULL) free(parameters->source_file);
        parameters->source_file = malloc(strlen(argv[i])*sizeof(char)+1);
        strcpy(parameters->source_file, argv[i]);
      }
      else print_error("Error reading command line input '-source_file'");
    }

    // Path to write histories to (-histories_file)
    else if(strcmp(arg, "-histories_file") == 0){
      if(++i < argc){
        if(parameters->histories_file != NULL) free(parameters->histories_file);
        parameters->histories_file = malloc(strlen(argv[i])*sizeof(char)+1);
        strcpy(parameters->histories_file, argv[i]);
      }
      else print_error("Error reading command line input '-histories_file'");
    }

    // Convergence method
    else if(strcmp(arg, "-ramp_up") == 0){
      if(++i < argc){
        if(strcasecmp(argv[i], "false") == 0)
          parameters->ramp_up = FALSE;
        else if(strcasecmp(argv[i], "true") == 0)
          parameters->ramp_up = TRUE;
        else
          print_error("Invalid convergence method: must be 'true' or 'false'");
      }
      else print_error("Error reading command line input '-ramp_up'");
    }

    // Unknown command line option
    else print_error("Error reading command line input");
  }

  // Validate Inputs
  if(parameters->entropy_bins <= 0)
    parameters->entropy_bins = ceil(pow(parameters->n_particles/20, 1.0/3.0));
  if(parameters->write_tally == TRUE && parameters->tally_file == NULL)
    parameters->tally_file = "tally.dat";
  if(parameters->write_entropy == TRUE && parameters->entropy_file == NULL)
    parameters->entropy_file = "entropy.dat";
  if(parameters->write_keff == TRUE && parameters->keff_file == NULL)
    parameters->keff_file = "keff.dat";
  if(parameters->write_bank == TRUE && parameters->bank_file == NULL)
    parameters->bank_file = "bank.dat";
  if(parameters->write_source == TRUE && parameters->source_file == NULL)
    parameters->source_file = "source.dat";
  if(parameters->write_histories == TRUE && parameters->histories_file == NULL)
    parameters->histories_file = "histories.dat";
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
  if(parameters->Lx <= 0 || parameters->Ly <= 0 || parameters->Lz <= 0)
    print_error("Length of domain must be positive in x, y, and z dimension");
  if(parameters->xs_f < 0 || parameters->xs_a < 0 || parameters->xs_s < 0)
    print_error("Macroscopic cross section values cannot be negative");

  return;
}

void read_convergence_parameters(Parameters *parameters)
{
  char line[512];
  char *s;
  FILE *fp;
  int i = 0;

  fp = fopen("convergence_parameters", "r");

  while((s = fgets(line, sizeof(line), fp)) != NULL){

    if(line[0] == '#' || line[0] == '\n') continue;

    // Number of stages
    s = strtok(line, " \n");
    parameters->cnvg_n_stages = atoi(s);

    break;
  }

  parameters->cnvg_n_particles = malloc(parameters->cnvg_n_stages*sizeof(int));
  parameters->cnvg_n_generations = malloc(parameters->cnvg_n_stages*sizeof(int));

  while((s = fgets(line, sizeof(line), fp)) != NULL){

    if(line[0] == '#' || line[0] == '\n') continue;

    // Number of particles in stage
    s = strtok(line, " ");
    parameters->cnvg_n_particles[i] = atoi(s);

    // Number of generations in stage
    s = strtok(NULL, " ");
    parameters->cnvg_n_generations[i] = atoi(s);

    i++;
  }

  fclose(fp);

  if(parameters->cnvg_n_particles[parameters->cnvg_n_stages-1] != parameters->n_particles){
    printf("WARNING: Number of particles in ramp-up does not match number of particles in parameters file.\n");
    parameters->n_particles = parameters->cnvg_n_particles[parameters->cnvg_n_stages-1];
  }
  if(parameters->n_active < parameters->n_batches){
    printf("WARNING: Setting number of active batches equal to number of total batches.\n"); 
    parameters->n_active = parameters->n_batches;
  }

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
  if(parameters->ramp_up == TRUE)
    printf("Convergence method:             Ramp up\n");
  else
    printf("Convergence method:             Default\n");
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

void init_output(Parameters *parameters)
{
  FILE *fp = NULL; // file pointer for output

  // Set up file to output tallies
  if(parameters->write_tally == TRUE){
    fp = fopen(parameters->tally_file, "w");
    fclose(fp);
  }

  // Set up file to output shannon entropy to assess source convergence
  if(parameters->write_entropy == TRUE){
    fp = fopen(parameters->entropy_file, "w");
    fclose(fp);
  }

  // Set up file to output keff
  if(parameters->write_keff == TRUE){
    fp = fopen(parameters->keff_file, "w");
    fclose(fp);
  }

  // Set up file to output particle bank
  if(parameters->write_bank == TRUE){
    fp = fopen(parameters->bank_file, "w");
    fclose(fp);
  }

  // Set up file to output source distribution
  if(parameters->write_source == TRUE){
    fp = fopen(parameters->source_file, "w");
    fclose(fp);
  }

  // Set up file to output total number of histories at each generation
  if(parameters->write_histories == TRUE){
    fp = fopen(parameters->source_file, "w");
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

void write_entropy(double H, char *filename)
{
  FILE *fp;

  fp = fopen(filename, "a");
  fprintf(fp, "%.10f\n", H);
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

void write_bank(Bank *b, char *filename)
{
  int i;
  FILE *fp;

  fp = fopen(filename, "a");

  for(i=0; i<b->n; i++){
    fprintf(fp, "%.10f %.10f %.10f ", b->p[i].x, b->p[i].y, b->p[i].z);
  }
  fprintf(fp, "\n");

  fclose(fp);
  return;
}

void write_histories(unsigned long n_histories, char *filename)
{
  FILE *fp;

  fp = fopen(filename, "a");
  fprintf(fp, "%lu\n", n_histories);
  fclose(fp);

  return;
}

void write_source(Parameters *parameters, Geometry *geometry, Bank *b, char *filename)
{
  int i, j, k;
  double dx, dy, dz;
  unsigned long ix, iy, iz;
  unsigned long l;
  unsigned long n;
  double *dist;
  Particle *p;
  FILE *fp;

  // Number of grid boxes in each dimension
  n = parameters->n_bins;

  // Find grid spacing
  dx = geometry->Lx/n;
  dy = geometry->Ly/n;
  dz = geometry->Lz/n;

  // Allocate array to keep track of number of sites in each grid box
  dist = calloc(n*n*n, sizeof(double));

  for(l=0; l<b->n; l++){
    p = &(b->p[l]);

    // Find the indices of the grid box of the particle
    ix = p->x/dx;
    iy = p->y/dy;
    iz = p->z/dz;

    dist[ix + n*iy + n*n*iz]++;
  }

  // Normalize by number of particles
  for(l=0; l<n*n*n; l++){
    dist[l] /= b->n;
  }

  fp = fopen(filename, "a");

  for(i=0; i<n; i++){
    for(j=0; j<n; j++){
      for(k=0; k<n; k++){
        fprintf(fp, "%e ", dist[i + n*j + n*n*k]);
      }
      fprintf(fp, "\n");
    }
  }

  fclose(fp);

  free(dist);

  return;
}

void load_source(Bank *b)
{
  unsigned long stat;
  FILE *fp;

  fp = fopen("source.dat", "rb");
  if(fp == NULL){
    print_error("Couldn't open source file.");
  }
  stat = fread(b->p, sizeof(Particle), b->sz, fp);
  if(stat != b->sz){
    print_error("Error loading source.");
  }
  fclose(fp);

  return;
}

void save_source(Bank *b)
{
  unsigned long stat;
  FILE *fp;

  fp = fopen("source.dat", "wb");
  stat = fwrite(b->p, sizeof(Particle), b->n, fp);
  if(stat != b->n){
    print_error("Error saving source.");
  }
  fclose(fp);

  return;
}
