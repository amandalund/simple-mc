#include "simple_mc.h"
#include "global.h"

Parameters *init_parameters(void)
{
  Parameters *p = malloc(sizeof(Parameters));

  p->n_particles = 1000000;
#ifdef _OPENMP
  p->n_threads = omp_get_num_procs();
#else
  p->n_threads = 1;
#endif
  p->n_batches = 10;
  p->n_generations = 1;
  p->n_active = 10;
  p->bc = REFLECT;
  p->n_nuclides = 1;
  p->tally = TRUE;
  p->n_bins = 16;
  p->seed = 1;
  p->nu = 2.5;
  p->xs_f = 0.012;
  p->xs_a = 0.03;
  p->xs_s = 0.27;
  p->Lx = 400;
  p->Ly = 400;
  p->Lz = 400;
  p->load_source = FALSE;
  p->save_source = FALSE;
  p->write_tally = FALSE;
  p->write_entropy = FALSE;
  p->write_msd = FALSE;
  p->write_keff = FALSE;
  p->write_bank = FALSE;
  p->write_source = FALSE;
  p->tally_file = NULL;
  p->entropy_file = NULL;
  p->msd_file = NULL;
  p->keff_file = NULL;
  p->bank_file = NULL;
  p->source_file = NULL;

  return p;
}

Geometry *init_geometry(Parameters *parameters)
{
  Geometry *g = malloc(sizeof(Geometry));

  g->Lx = parameters->Lx;
  g->Ly = parameters->Ly;
  g->Lz = parameters->Lz;
  g->bc = parameters->bc;

  return g;
}

Tally *init_tally(Parameters *parameters)
{
  Tally *t = malloc(sizeof(Tally));

  t->tallies_on = FALSE;
  t->n = parameters->n_bins;
  t->dx = parameters->Lx/t->n;
  t->dy = parameters->Ly/t->n;
  t->dz = parameters->Lz/t->n;
  t->flux = calloc(t->n*t->n*t->n, sizeof(double));

  return t;
}

Material *init_material(Parameters *parameters)
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

Bank *init_source_bank(Parameters *parameters, Geometry *geometry)
{
  unsigned long i_p; // index over particles
  Bank *source_bank;

  // Initialize source bank
  source_bank = init_bank(parameters->n_particles);

  // Sample source particles or load a source
  if(parameters->load_source == TRUE){
    load_source(source_bank);
    source_bank->n = parameters->n_particles;
  }
  else{
    for(i_p=0; i_p<parameters->n_particles; i_p++){
      sample_source_particle(geometry, &(source_bank->p[i_p]));
      source_bank->n++;
    }
  }

  return source_bank;
}

void init_fission_bank(Parameters *parameters)
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

void sample_source_particle(Geometry *geometry, Particle *p)
{
  p->alive = TRUE;
  p->energy = 1;
  p->last_energy = 1;
  p->mu = rn()*2 - 1;
  p->phi = rn()*2*PI;
  p->u = p->mu;
  p->v = sqrt(1 - p->mu*p->mu)*cos(p->phi);
  p->w = sqrt(1 - p->mu*p->mu)*sin(p->phi);
  p->x = rn()*geometry->Lx;
  p->y = rn()*geometry->Ly;
  p->z = rn()*geometry->Lz;

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
