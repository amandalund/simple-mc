#ifndef HEADER
#define HEADER

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<sys/time.h>
#include<math.h>
#include<float.h>
#include<unistd.h>
#include<string.h>

typedef struct Parameters_{
  unsigned long long seed; // RNG seed
  unsigned long n_particles; // number of particles
  int n_batches; // number of batches
  int n_generations; // number of generations per batch
  int n_active; // number of active batches
  int bc; // boundary conditions
  int n_nuclides; // number of nuclides in material
  int tally; // whether to tally
  int n_bins; // number of bins in each dimension of mesh
  double nu; // average number of fission neutrons produced
  double xs_a; // absorption macro xs
  double xs_s; // scattering macro xs
  double xs_f; // fission macro xs
  double gx; // geometry size in x
  double gy; // geometry size in y
  double gz; // geometry size in z
  int load_source; // load the source bank from source.dat
  int save_source; // save the source bank at end of simulation
  int write_tally; // whether to output tallies
  int write_entropy; // whether to output shannon entropy
  int write_keff; // whether to output keff
  int write_bank; // whether to output particle bank
  int write_source; // whether to output source distribution
  char *tally_file; // path to write tallies to
  char *entropy_file; // path to write shannon entropy to
  char *keff_file; // path to write keff to
  char *bank_file; // path to write particle bank to
  char *source_file; // path to write source distribution to
} Parameters;

typedef struct Particle_{
  int alive;
  double energy;
  double last_energy;
  double mu; // cosine of polar angle
  double phi; // azimuthal angle
  double u; // direction
  double v;
  double w;
  double x; // position
  double y;
  double z;
  int surface_crossed;
  int event;
} Particle;

typedef struct Geometry_{
  int bc;
  double x;
  double y;
  double z;
} Geometry;

typedef struct Nuclide_{
  double xs_f; // fission micro xs
  double xs_a; // absorption micro xs
  double xs_s; // scattering micro xs
  double xs_t; // total micro xs
  double atom_density; // atomic density of nuclide in material
} Nuclide;

typedef struct Material_{
  double xs_f; // fission macro xs
  double xs_a; // absorption macro xs
  double xs_s; // scattering macro xs
  double xs_t; // total macro xs
  int n_nuclides;
  Nuclide *nuclides;
} Material;

typedef struct Tally_{
  int tallies_on; // whether tallying is currently turned on
  int n; // mumber of grid boxes in each dimension 
  double dx; // grid spacing
  double dy;
  double dz;
  double *flux;
} Tally;

typedef struct Bank_{
  unsigned long n; // number of particles
  unsigned long sz; // size of bank
  Particle *p; // particle array
  void (*resize)(struct Bank_ *b);
} Bank;

// io.c function prototypes
void parse_parameters(void);
void read_CLI(int argc, char *argv[]);
void print_error(char *message);
void print_parameters(void);
void border_print(void);
void fancy_int(long a);
void center_print(const char *s, int width);
void init_output(void);
void write_tally(Tally *t, char *filename);
void write_entropy(double H, char *filename);
void write_keff(double *keff, int n, char *filename);
void write_bank(Bank *b, char *filename);
void write_source(Bank *b, char *filename);
void load_source(Bank *b);
void save_source(Bank *b);

// utils.c funtion prototypes
double timer(void);

// prng.c function prototypes
double rn(void);
int rni(int min, int max);
void set_stream(int rn_stream);
void set_initial_seed(unsigned long long rn_seed0);
void rn_skip(long long n);

// initialize.c function prototypes
void init_problem(int argc, char *argv[]);
Parameters *init_parameters(void);
Geometry *init_geometry(void);
Tally *init_tally(void);
Material *init_material(void);
void init_fission_bank(void);
void init_source_bank(void);
Bank *init_bank(unsigned long n_particles);
void sample_source_particle(Particle *p);
void sample_fission_particle(Particle *p, Particle *p_old);
void copy_particle(Particle *dest, Particle *source);
void resize_particles(Bank *b);
void free_bank(Bank *b);
void free_material(Material *m);
void free_tally(Tally *t);
void free_problem(void);

// transport.c function prototypes
void transport(Particle *p);
void calculate_xs(void);
double distance_to_boundary(Particle *p);
double distance_to_collision(void);
void cross_surface(Particle *p);
void collision(Particle *p);

// eigenvalue.c function prototypes
void run_eigenvalue(void);
void synchronize_bank(void);
double shannon_entropy(Bank *b);
void calculate_keff(double *mean, double *std, int n);

// tally.c function prototypes
void score_tally(Tally *t, Particle *p);

#endif
