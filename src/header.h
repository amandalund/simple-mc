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
#include<omp.h>

typedef struct Parameters_{
  unsigned long long seed;   // RNG seed
  unsigned long n_particles; // number of particles
  int n_threads;             // number of openmp threads
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
void parse_params();
void read_CLI(int argc, char *argv[]);
void print_error(char *message);
void print_params();
void border_print(void);
void fancy_int(long a);
void center_print(const char *s, int width);
void init_output(Parameters *params, FILE *fp);
void write_tally(Tally *t, FILE *fp, char *filename);
void write_entropy(double H, FILE *fp, char *filename);
void write_keff(double *keff, int n, FILE *fp, char *filename);
void write_bank(Bank *b, FILE *fp, char *filename);
void write_source(Geometry *g, Bank *b, Parameters *params, FILE *fp, char *filename);
void load_source(Bank *b);
void save_source(Bank *b);

// utils.c funtion prototypes
double timer(void);

// prng.c function prototypes
double rn(unsigned long long *seed);
int rni(unsigned long long *seed, int min, int max);
unsigned long long rn_skip(unsigned long long seed, long long n);

// initialize.c function prototypes
void init_problem(int argc, char *argv[]);
Parameters *init_params(void);
Geometry *init_geometry(Parameters *params);
Tally *init_tally(Parameters *params);
Material *init_material(Parameters *params, unsigned long long *seed);
Bank *init_bank(unsigned long n_particles);
void sample_source_particle(Particle *p, Geometry *g, unsigned long long *seed);
void sample_fission_particle(Particle *p, Particle *p_old, unsigned long long *seed);
void copy_particle(Particle *dest, Particle *source);
void resize_particles(Bank *b);
void free_bank(Bank *b);
void free_material(Material *m);
void free_tally(Tally *t);

// transport.c function prototypes
void transport(Particle *p, Geometry *g, Material *m, Tally *t, Parameters *params, unsigned long long *seed);
void calculate_xs(Particle *p, Material *m);
double distance_to_boundary(Particle *p, Geometry *g);
double distance_to_collision(Material *m, unsigned long long *seed);
void cross_surface(Particle *p, Geometry *g);
void collision(Particle *p, Material *m, double nu, unsigned long long *seed);

// eigenvalue.c function prototypes
void merge_fission_banks();
void synchronize_bank(Geometry *g, unsigned long long *seed);
double shannon_entropy(Geometry *g, Bank *b);
void calculate_keff(double *keff, double *mean, double *std, int n);

// tally.c function prototypes
void score_tally(Tally *t, Particle *p, Material *m, Parameters *params);

#endif
