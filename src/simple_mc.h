#ifndef SIMPLE_MC
#define SIMPLE_MC

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
  unsigned long long seed; // RNG seed
  unsigned long n_particles; // number of particles
  unsigned long n_histories; // cumulative number of histories
  int n_threads; // number of openmp threads
  int n_batches; // number of batches
  int n_generations; // number of generations per batch
  int n_active; // number of active batches
  int bc; // boundary conditions
  int source; // initial source distribution
  int n_nuclides; // number of nuclides in material
  int tally; // whether to tally
  int n_bins; // number of bins in each dimension of mesh
  double nu; // average number of fission neutrons produced
  double xs_a; // absorption macro xs
  double xs_s; // scattering macro xs
  double xs_f; // fission macro xs
  double Lx; // domain length in x
  double Ly; // domain length in y
  double Lz; // domain length in z
  int load_source; // load the source bank from source.dat
  int save_source; // save the source bank at end of simulation
  int write_tally; // whether to output tallies
  int write_entropy; // whether to output shannon entropy
  int write_histories; // whether to output histories
  int write_msd; // whether to output mean-squared distance
  int write_keff; // whether to output keff
  int write_bank; // whether to output particle bank
  int write_source; // whether to output source distribution
  char *tally_file; // path to write tallies to
  char *entropy_file; // path to write shannon entropy to
  char *histories_file; // path to write histories to
  char *msd_file; // path to write mean-squared distance to
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
  double Lx;
  double Ly;
  double Lz;
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
void parse_parameters(Parameters *parameters);
void read_CLI(int argc, char *argv[], Parameters *parameters);
void print_error(char *message);
void print_parameters(Parameters *parameters);
void border_print(void);
void fancy_int(long a);
void center_print(const char *s, int width);
void init_output(Parameters *parameters);
void write_tally(Tally *t, char *filename);
void write_entropy(double H, char *filename);
void write_histories(unsigned long histories, char *filename);
void write_msd(double msd, char *filename);
void write_keff(double *keff, int n, char *filename);
void write_bank(Bank *b, char *filename);
void write_source(Parameters *parameters, Geometry *geometry, Bank *b, char *filename);
void load_source(Bank *b);
void save_source(Bank *b);

// utils.c funtion prototypes
double timer(void);
void copy_particle(Particle *dest, Particle *source);

// prng.c function prototypes
double rn(void);
int rni(int min, int max);
void set_stream(int rn_stream);
void set_initial_seed(unsigned long long rn_seed0);
void rn_skip(long long n);

// initialize.c function prototypes
Parameters *init_parameters(void);
Geometry *init_geometry(Parameters *parameters);
Tally *init_tally(Parameters *parameters);
Material *init_material(Parameters *parameters);
void init_fission_bank(Parameters *parameters);
Bank *init_source_bank(Parameters *parameters, Geometry *geometry);
Bank *init_bank(unsigned long n_particles);
void sample_source_particle(Geometry *geometry, Particle *p);
void resize_particles(Bank *b);
void free_bank(Bank *b);
void free_material(Material *m);
void free_tally(Tally *t);

// transport.c function prototypes
void transport(Parameters *parameters, Geometry *geometry, Material *material, Bank *source_bank, Tally *tally, Particle *p);
void calculate_xs(Material *material);
double distance_to_boundary(Geometry *geometry, Particle *p);
double distance_to_collision(Material *material);
void cross_surface(Geometry *geometry, Particle *p);
void collision(Material *material, double nu, Particle *p);
void sample_fission_particle(Particle *p, Particle *p_old);

// eigenvalue.c function prototypes
void run_eigenvalue(Parameters *parameters, Geometry *geometry, Material *material, Bank *source_bank, Tally *tally, double *keff);
void merge_fission_banks(void);
void synchronize_bank(Bank *source_bank);
double shannon_entropy(Geometry *geometry, Bank *b);
double mean_squared_distance(Bank *b);
void calculate_keff(double *keff, double *mean, double *std, int n);

// tally.c function prototypes
void score_tally(Parameters *parameters, Material *material, Tally *t, Particle *p);

#endif
