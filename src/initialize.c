#include "header.h"
#include "global.h"

void init_problem(int argc, char *argv[])
{
  // Get inputs: set parameters to default values, then parse parameter file,
  // then override with any command line inputs
  parameters = init_parameters();
  parse_parameters();
  read_CLI(argc, argv);
  print_parameters();

  // Set initial RNG seed
  set_initial_seed(parameters->seed);
  set_stream(STREAM_INIT);

  // Set OpenMP specific variables
#ifdef _OPENMP
  omp_set_num_threads(parameters->n_threads);
#pragma omp parallel
{
  n_threads = omp_get_num_threads();
  thread_id = omp_get_thread_num();
}
#endif

  // Set up output files
  init_output();

  // Set up geometry
  geometry = init_geometry();

  // Set up material
  material = init_material();

  // Set up tallies
  tally = init_tally();

  // Set up fission banks
  init_fission_bank();

  // Set up source bank and initial source distribution
  init_source_bank();

  // Set up array for keff
  keff = calloc(parameters->n_active, sizeof(double));

  return;
}

Parameters *init_parameters()
{
  Parameters *parameters = malloc(sizeof(Parameters));

#ifdef _OPENMP
  parameters->n_threads = omp_get_num_procs();
#else
  parameters->n_threads = 1;
#endif
  parameters->n_particles = 1000000;
  parameters->n_batches = 10;
  parameters->n_generations = 1;
  parameters->n_active = 10;
  parameters->bc = REFLECT;
  parameters->n_nuclides = 1;
  parameters->tally = TRUE;
  parameters->n_bins = 16;
  parameters->seed = 1;
  parameters->nu = 2.5;
  parameters->xs_f = 0.012;
  parameters->xs_a = 0.03;
  parameters->xs_s = 0.27;
  parameters->gx = 400;
  parameters->gy = 400;
  parameters->gz = 400;
  parameters->load_source = FALSE;
  parameters->save_source = FALSE;
  parameters->write_tally = FALSE;
  parameters->write_entropy = FALSE;
  parameters->write_keff = FALSE;
  parameters->write_bank = FALSE;
  parameters->write_source = FALSE;
  parameters->tally_file = NULL;
  parameters->entropy_file = NULL;
  parameters->keff_file = NULL;
  parameters->bank_file = NULL;
  parameters->source_file = NULL;

  return parameters;
}

Geometry *init_geometry(void)
{
  Geometry *g = malloc(sizeof(Geometry));

  g->x = parameters->gx;
  g->y = parameters->gy;
  g->z = parameters->gz;
  g->bc = parameters->bc;

  return g;
}

Tally *init_tally(void)
{
  Tally *t = malloc(sizeof(Tally));

  t->tallies_on = FALSE;
  t->n = parameters->n_bins;
  t->dx = parameters->gx/t->n;
  t->dy = parameters->gy/t->n;
  t->dz = parameters->gz/t->n;
  t->flux = calloc(t->n*t->n*t->n, sizeof(double));

  return t;
}

Material *init_material(void)
{
  int i;
  Nuclide sum = {0, 0, 0, 0, 0};

  // Hardwire the material macroscopic cross sections for now to produce a keff
  // close to 1 (fission, absorption, scattering, total, atomic density)
  Nuclide macro = {parameters->xs_f, parameters->xs_a, parameters->xs_s,
     parameters->xs_f + parameters->xs_a + parameters->xs_s, 1.0};

  Material *m = malloc(sizeof(Material));
  m->n_nuclides = parameters->n_nuclides;
  m->nuclides = malloc(m->n_nuclides*sizeof(Nuclide));

  // Generate some arbitrary microscopic cross section values and atomic
  // densities for each nuclide in the material such that the total macroscopic
  // cross sections evaluate to what is hardwired above
  for(i=0; i<m->n_nuclides; i++){
    if(i<m->n_nuclides-1){
      m->nuclides[i].atom_density = rn()*macro.atom_density;
      macro.atom_density -= m->nuclides[i].atom_density;
    }
    else{
      m->nuclides[i].atom_density = macro.atom_density;
    }
    m->nuclides[i].xs_a = rn();
    sum.xs_a += m->nuclides[i].xs_a * m->nuclides[i].atom_density;
    m->nuclides[i].xs_f = rn();
    sum.xs_f += m->nuclides[i].xs_f * m->nuclides[i].atom_density;
    m->nuclides[i].xs_s = rn();
    sum.xs_s += m->nuclides[i].xs_s * m->nuclides[i].atom_density;
  }
  for(i=0; i<m->n_nuclides; i++){
    m->nuclides[i].xs_a /= sum.xs_a/macro.xs_a;
    m->nuclides[i].xs_f /= sum.xs_f/macro.xs_f;
    m->nuclides[i].xs_s /= sum.xs_s/macro.xs_s;
    m->nuclides[i].xs_t = m->nuclides[i].xs_a + m->nuclides[i].xs_s;
  }

  m->xs_f = parameters->xs_f;
  m->xs_a = parameters->xs_a;
  m->xs_s = parameters->xs_s;
  m->xs_t = parameters->xs_a + parameters->xs_s;

  return m;
}

void init_source_bank(void)
{
  unsigned long i_p; // index over particles

  // Initialize source bank
  source_bank = init_bank(parameters->n_particles);

  // Sample source particles or load a source
  if(parameters->load_source == TRUE){
    load_source(source_bank);
    source_bank->n = parameters->n_particles;
  }
  else{
    for(i_p=0; i_p<parameters->n_particles; i_p++){
      sample_source_particle(&(source_bank->p[i_p]));
      source_bank->n++;
    }
  }

  return;
}

void init_fission_bank(void)
{
  // Allocate one fission bank for each thread and one master fission bank to
  // collect fission sites at end of each generation
#ifdef _OPENMP
#pragma omp parallel
{
  if(thread_id == 0){
    fission_bank = init_bank(2*parameters->n_particles);
  }
  else{
    fission_bank = init_bank(2*parameters->n_particles/n_threads);
  }
}
  master_fission_bank = init_bank(2*parameters->n_particles);
#else
  fission_bank = init_bank(2*parameters->n_particles);
#endif

  return;
}

Bank *init_bank(unsigned long n_particles)
{
  Bank *b = malloc(sizeof(Bank));
  b->p = malloc(n_particles*sizeof(Particle));
  b->sz = n_particles;
  b->n = 0;
  b->resize = resize_particles;

  return b;
}

void sample_source_particle(Particle *p)
{
  p->alive = TRUE;
  p->energy = 1;
  p->last_energy = 1;
  p->mu = rn()*2 - 1;
  p->phi = rn()*2*PI;
  p->u = p->mu;
  p->v = sqrt(1 - p->mu*p->mu)*cos(p->phi);
  p->w = sqrt(1 - p->mu*p->mu)*sin(p->phi);
  p->x = rn()*geometry->x;
  p->y = rn()*geometry->y;
  p->z = rn()*geometry->z;

  return;
}

void sample_fission_particle(Particle *p, Particle *p_old)
{
  p->alive = TRUE;
  p->energy = 1;
  p->last_energy = 1;
  p->mu = rn()*2 - 1;
  p->phi = rn()*2*PI;
  p->u = p->mu;
  p->v = sqrt(1 - p->mu*p->mu)*cos(p->phi);
  p->w = sqrt(1 - p->mu*p->mu)*sin(p->phi);
  p->x = p_old->x;
  p->y = p_old->y;
  p->z = p_old->z;

  return;
}

void copy_particle(Particle *dest, Particle *source)
{
  dest->alive = source->alive;
  dest->energy = source->energy;
  dest->last_energy = source->last_energy;
  dest->mu = source->mu;
  dest->phi = source->phi;
  dest->u = source->u;
  dest->v = source->v;
  dest->w = source->w;
  dest->x = source->x;
  dest->y = source->y;
  dest->z = source->z;
  dest->event = source->event;

  return;
}

void resize_particles(Bank *b)
{
  b->p = realloc(b->p, sizeof(Particle)*2*b->sz);
  b->sz = 2*b->sz;

  return;
}

void free_bank(Bank *b)
{
  free(b->p);
  b->p = NULL;
  free(b);
  b = NULL;

  return;
}

void free_material(Material *m)
{
  free(m->nuclides);
  m->nuclides = NULL;
  free(m);
  m = NULL;

  return;
}

void free_tally(Tally *t)
{
  free(t->flux);
  t->flux = NULL;
  free(t);
  t = NULL;

  return;
}

void free_problem(void)
{
#ifdef _OPENMP
#pragma omp parallel
{
  free_bank(fission_bank);
}
  free_bank(master_fission_bank);
#else
  free_bank(fission_bank);
#endif
  free_bank(source_bank);
  free_tally(tally);
  free_material(material);
  free(geometry);
  free(keff);
  free(parameters);

  return;
}
