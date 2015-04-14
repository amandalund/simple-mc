#include "simple_transport_header.h"

Input *set_default_input(void)
{
  Input *I = malloc(sizeof(Input));

  I->n_particles = 10000;
  I->n_batches = 20;
  I->n_generations = 1;
  I->n_active = 10;
  I->bc = REFLECT;
  I->n_nuclides = 60;

  return I;
}

Geometry *init_geometry(int bc)
{
  Geometry *g = malloc(sizeof(Geometry));

  g->x = 1000;
  g->y = 1000;
  g->z = 1000;
  g->bc = bc;
  g->surface_crossed = -1;

  return g;
}

Material *init_material(int n_nuclides)
{
  int i;
  Nuclide sum = {0, 0, 0, 0, 0};

  // Hardwire the material macroscopic cross sections for now to produce a keff
  // close to 1 (fission, absorption, elastic, total, atomic density)
  Nuclide macro = {2.29, 3.42, 2.29, 8.00, 1.0};

  Material *m = malloc(sizeof(Material));
  m->n_nuclides = n_nuclides;
  m->nuclides = malloc(n_nuclides*sizeof(Nuclide));

  // Generate some arbitrary microscopic cross section values and atomic
  // densities for each nuclide in the material such that the total macroscopic
  // cross sections evaluate to what is hardwired above
  for(i=0; i<n_nuclides; i++){
    if(i<n_nuclides-1){
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
    m->nuclides[i].xs_e = rn();
    sum.xs_e += m->nuclides[i].xs_e * m->nuclides[i].atom_density;
  }
  for(i=0; i<n_nuclides; i++){
    m->nuclides[i].xs_a /= sum.xs_a/macro.xs_a;
    m->nuclides[i].xs_f /= sum.xs_f/macro.xs_f;
    m->nuclides[i].xs_e /= sum.xs_e/macro.xs_e;
    m->nuclides[i].xs_t = m->nuclides[i].xs_a + m->nuclides[i].xs_f + m->nuclides[i].xs_e;
  }

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
