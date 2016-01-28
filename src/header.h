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

#define TRUE 1
#define FALSE 0

// Constants
#define PI 3.1415926535898
#define D_INF DBL_MAX

// Geometry boundary conditions
#define VACUUM 0
#define REFLECT 1
#define PERIODIC 2

// Reaction types
#define TOTAL 0
#define ABSORPTION 1
#define SCATTER 2
#define FISSION 3

// Surfaces
#define X0 0
#define X1 1
#define Y0 2
#define Y1 3
#define Z0 4
#define Z1 5

typedef struct Parameters_{
  unsigned long seed;        // RNG seed
  unsigned long n_particles; // number of particles
  int n_batches;             // number of batches
  int n_generations;         // number of generations per batch
  int n_active;              // number of active batches
  int bc;                    // boundary conditions
  int n_nuclides;            // number of nuclides in material
  int tally;                 // whether to tally
  int n_bins;                // number of bins in each dimension of mesh
  double nu;                 // average number of fission neutrons produced
  double xs_a;               // absorption macro xs
  double xs_s;               // scattering macro xs
  double xs_f;               // fission macro xs
  double gx;                 // geometry size in x
  double gy;                 // geometry size in y
  double gz;                 // geometry size in z
  int load_source;           // load the source bank from source.dat
  int save_source;           // save the source bank at end of simulation
  int write_tally;           // whether to output tallies
  int write_entropy;         // whether to output shannon entropy
  int write_keff;            // whether to output keff
  int write_bank;            // whether to output particle bank
  int write_source;          // whether to output source distribution
  char *tally_file;          // path to write tallies to
  char *entropy_file;        // path to write shannon entropy to
  char *keff_file;           // path to write keff to
  char *bank_file;           // path to write particle bank to
  char *source_file;         // path to write source distribution to
} Parameters;

typedef struct Particle_{
  int alive;
  double energy;
  double last_energy;
  double mu;                 // cosine of polar angle
  double phi;                // azimuthal angle
  double u;                  // direction
  double v;
  double w;
  double x;                  // position
  double y;
  double z;
  int event;
} Particle;

typedef struct Geometry_{
  int bc;
  double x;
  double y;
  double z;
  int surface_crossed;
} Geometry;

typedef struct Nuclide_{
  double xs_f;               // fission micro xs
  double xs_a;               // absorption micro xs
  double xs_s;               // scattering micro xs
  double xs_t;               // total micro xs
  double atom_density;       // atomic density of nuclide in material
} Nuclide;

typedef struct Material_{
  double xs_f;               // fission macro xs
  double xs_a;               // absorption macro xs
  double xs_s;               // scattering macro xs
  double xs_t;               // total macro xs
  int n_nuclides;
  Nuclide *nuclides;
} Material;

typedef struct Tally_{
  int tallies_on;            // whether tallying is currently turned on
  int n;                     // mumber of grid boxes in each dimension 
  double dx;                 // grid spacing
  double dy;
  double dz;
  double *flux;
} Tally;

typedef struct Bank_{
  unsigned long n;           // number of particles
  unsigned long sz;          // size of bank
  Particle *p;               // particle array
  void (*resize)(struct Bank_ *b);
} Bank;

// io.c function prototypes
void parse_params(char *filename, Parameters *params);
void read_CLI(int argc, char *argv[], Parameters *params);
void print_error(char *message);
void print_params(Parameters *params);
void border_print(void);
void fancy_int(long a);
void center_print(const char *s, int width);
void write_tally(Tally *t, FILE *fp, char *filename);
void write_entropy(double H, FILE *fp, char *filename);
void write_keff(double *keff, int n, FILE *fp, char *filename);
void write_bank(Bank *b, FILE *fp, char *filename);
void write_source(Geometry *g, Bank *b, Parameters *params, FILE *fp, char *filename);
void load_source(Bank *b);
void save_source(Bank *b);

// utils.c funtion prototypes
double rn(unsigned long *seed);
int rni(unsigned long *seed, int min, int max);
double timer(void);

// initialize.c function prototypes
Parameters *set_default_params(void);
void init_output(Parameters *params, FILE *fp);
Geometry *init_geometry(Parameters *params);
Tally *init_tally(Parameters *params);
Material *init_material(Parameters *params);
Bank *init_bank(unsigned long n_particles);
void sample_source_particle(Particle *p, Geometry *g, Parameters *params);
void sample_fission_particle(Particle *p, Particle *p_old, Parameters *params);
void resize_particles(Bank *b);
void free_bank(Bank *b);
void free_material(Material *m);
void free_tally(Tally *t);

// transport.c function prototypes
void transport(Particle *p, Geometry *g, Material *m, Tally *t, Bank *fission_bank, double keff, Parameters *params);
void calculate_xs(Particle *p, Material *m);
double distance_to_boundary(Particle *p, Geometry *g);
double distance_to_collision(Material *m, Parameters *params);
void cross_surface(Particle *p, Geometry *g);
void collision(Particle *p, Material *m, Bank *fission_bank, double keff, Parameters *params);

// eigenvalue.c function prototypes
void synchronize_bank(Bank *source_bank, Bank *fission_bank, Geometry *g, Parameters *params);
double shannon_entropy(Geometry *g, Bank *b, Parameters *params);
void calculate_keff(double *keff, double *mean, double *std, int n);

// tally.c function prototypes
void score_tally(Tally *t, Particle *p, Material *m, Parameters *params);

#endif
