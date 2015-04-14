#ifndef SIMPLE_TRANSPORT_HEADER
#define SIMPLE_TRANSPORT_HEADER

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
#define PERIODIC 3

// Surfaces
#define X 0
#define Y 1
#define Z 2

// Inputs
typedef struct Input_{
  unsigned long n_particles; // number of particles
  int n_batches;             // number of batches
  int n_generations;         // number of generations per batch
  int n_active;              // number of active batches
  int bc;                    // boundary conditions
  int n_nuclides;            // number of nuclides in material
} Input;

// Particle
typedef struct Particle_{
  int alive;
  double energy;
  double last_energy;
  double mu;      // cosine of polar angle
  double phi;     // azimuthal angle
  double u;       // direction
  double v;
  double w;
  double x;       // position
  double y;
  double z;
} Particle;

// Box geometry
typedef struct Geometry_{
  int bc;
  double x;
  double y;
  double z;
  int surface_crossed;
} Geometry;

typedef struct Nuclide_{
  double xs_f;         // Fission micro xs
  double xs_a;         // Absorption micro xs
  double xs_e;         // Elastic micro xs
  double xs_t;         // Total micro xs
  double atom_density; // Atomic density of nuclide in material
} Nuclide;

// Material
typedef struct Material_{
  double xs_f; // Fission macro xs
  double xs_a; // Absorption macro xs
  double xs_e; // Elastic macro xs
  double xs_t; // Total macro xs
  int n_nuclides;
  Nuclide *nuclides;
} Material;

// Particle bank
typedef struct Bank_{
  unsigned long n;  // number of particles
  unsigned long sz; // size of bank
  Particle *p;      // particle array
  void (*resize)(struct Bank_ *b);
} Bank;

// io.c function prototypes
void read_CLI(int argc, char *argv[], Input *input);
void print_CLI_error(void);
void print_inputs(Input *input);
void border_print(void);
void fancy_int(long a);
void center_print(const char *s, int width);
void print_particle(Particle *p);
void print_bank(Bank *b);

// utils.c funtion prototypes
//double rn(unsigned long * seed);
double rn(void);
double timer(void);

// initialize.c function prototypes
Input *set_default_input(void);
Geometry *init_geometry(int bc);
Material *init_material(int n_nuclides);
Bank *init_bank(unsigned long n_particles);
void sample_source_particle(Particle *p, Geometry *g);
void resize_particles(Bank *b);
void free_bank(Bank *b);
void free_material(Material *m);

// transport.c function prototypes
void transport(Particle *p, Geometry *g, Material *m, Bank *fission_bank);
void calculate_xs(Particle *p, Material *m);
double distance_to_boundary(Particle *p, Geometry *g);
double distance_to_collision(Material *m);
void cross_surface(Particle *p, Geometry *g);
void collision(Particle *p, Material *m, Bank *fission_bank);

// eigenvalue.c function prototypes
void synchronize_bank(Bank *source_bank, Bank *fission_bank, Geometry *g);
double shannon_entropy(Geometry *g, Bank *b);

#endif
