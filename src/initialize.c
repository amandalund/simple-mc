#include "header.h"

Parameters *set_default_params(void)
{
  Parameters *params = malloc(sizeof(Parameters));

  params->n_particles = 10000;
  params->n_batches = 20;
  params->n_generations = 1;
  params->n_active = 10;
  params->bc = REFLECT;
  params->n_nuclides = 60;
  params->tally = FALSE;
  params->n_bins = 10;
  params->seed = 1;
  params->xs_f = 2.29;
  params->xs_a = 3.42;
  params->xs_s = 2.29;
  params->gx = 1000;
  params->gy = 1000;
  params->gz = 1000;
  params->load_source = FALSE;
  params->save_source = FALSE;
  params->write_tally = FALSE;
  params->write_entropy = FALSE;
  params->write_keff = FALSE;
  params->write_bank = FALSE;
  params->tally_file = NULL;
  params->entropy_file = NULL;
  params->keff_file = NULL;
  params->bank_file = NULL;

  return params;
}

void init_output(Parameters *params, FILE *fp)
{
  // Set up file to output tallies
  if(params->write_tally == TRUE){
    fp = fopen(params->tally_file, "w");
    fclose(fp);
  }

  // Set up file to output shannon entropy to assess source convergence
  if(params->write_entropy == TRUE){
    fp = fopen(params->entropy_file, "w");
    fclose(fp);
  }

  // Set up file to output keff
  if(params->write_keff == TRUE){
    fp = fopen(params->keff_file, "w");
    fclose(fp);
  }

  // Set up file to output particle bank
  if(params->write_bank == TRUE){
    fp = fopen(params->bank_file, "w");
    fclose(fp);
  }

  return;
}

Geometry *init_geometry(Parameters *params)
{
  Geometry *g = malloc(sizeof(Geometry));

  g->x = params->gx;
  g->y = params->gy;
  g->z = params->gz;
  g->bc = params->bc;
  g->surface_crossed = -1;

  return g;
}

Tally *init_tally(Parameters *params)
{
  Tally *t = malloc(sizeof(Tally));

  t->tallies_on = FALSE;
  t->n = params->n_bins;
  t->dx = params->gx/t->n;
  t->dy = params->gy/t->n;
  t->dz = params->gz/t->n;
  t->flux = calloc(t->n*t->n*t->n, sizeof(double));

  return t;
}

Material *init_material(Parameters *params)
{
  int i;
  Nuclide sum = {0, 0, 0, 0, 0};

  // Hardwire the material macroscopic cross sections for now to produce a keff
  // close to 1 (fission, absorption, scattering, total, atomic density)
  Nuclide macro = {params->xs_f, params->xs_a, params->xs_s,
     params->xs_f + params->xs_a + params->xs_s, 1.0};

  Material *m = malloc(sizeof(Material));
  m->n_nuclides = params->n_nuclides;
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
    m->nuclides[i].xs_t = m->nuclides[i].xs_a + m->nuclides[i].xs_f + m->nuclides[i].xs_s;
  }

  m->xs_f = params->xs_f;
  m->xs_a = params->xs_a;
  m->xs_s = params->xs_s;
  m->xs_t = params->xs_f + params->xs_a + params->xs_f;

  return m;
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

void sample_source_particle(Particle *p, Geometry *g)
{
  p->alive = TRUE;
  p->energy = 1;
  p->last_energy = 0;
  p->mu = rn()*2 - 1;
  p->phi = rn()*2*PI;
  p->u = p->mu;
  p->v = sqrt(1 - p->mu*p->mu)*cos(p->phi);
  p->w = sqrt(1 - p->mu*p->mu)*sin(p->phi);
  p->x = rn()*g->x;
  p->y = rn()*g->y;
  p->z = rn()*g->z;

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
