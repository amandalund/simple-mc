#include "simple_transport_header.h"

// Main logic to move particle
void transport(Particle *p, Geometry *g, Material *m, Tally *t, Bank *fission_bank)
{
  while(p->alive){

    // Recalculate macro xs if particle has changed energy
    if(p->energy != p->last_energy){
      calculate_xs(p, m);
    }

    // Find distance to boundary
    double d_b = distance_to_boundary(p, g);

    // Find distance to collision
    double d_c = distance_to_collision(m);

    // Take smaller of two distances
    double d = d_b < d_c ? d_b : d_c;

    // Advance particle
    p->x = p->x + d*p->u;
    p->y = p->y + d*p->v;
    p->z = p->z + d*p->w;

    // Case where particle crosses boundary
    if(d_b < d_c){
      cross_surface(p, g);
    }
    // Case where particle has collision
    else{
      collision(p, m, fission_bank);

      // Score tallies
      if(t->tallies_on == TRUE){
        score_tally(t, p);
      }
    }
  }
  return;
}

// Calculates the macroscopic cross section of the material the particle is
// traveling through
void calculate_xs(Particle *p, Material *m)
{
  int i;

  // Reset macroscopic cross sections to 0
  m->xs_t = 0.0;
  m->xs_f = 0.0;
  m->xs_a = 0.0;
  m->xs_e = 0.0;

  for(i=0; i<m->n_nuclides; i++){

    Nuclide nuc = m->nuclides[i];

    // Add contribution from this nuclide to total macro xs
    m->xs_t += nuc.atom_density * nuc.xs_t;

    // Add contribution from this nuclide to fission macro xs
    m->xs_f += nuc.atom_density * nuc.xs_f;

    // Add contribution from this nuclide to absorption macro xs
    m->xs_a += nuc.atom_density * nuc.xs_a;

    // Add contribution from this nuclide to elastic macro xs
    m->xs_e += nuc.atom_density * nuc.xs_e;
  }

  return;
}

// Returns the distance to the nearest boundary for a particle traveling in a
// certain direction
double distance_to_boundary(Particle *p, Geometry *g)
{
  int i;
  double dist;
  double d = D_INF;
  int    surfaces[6] = {X, X, Y, Y, Z, Z};
  double p_angles[6] = {p->u, p->u, p->v, p->v, p->w, p->w};
  double p_coords[6] = {p->x, p->x, p->y, p->y, p->z, p->z};
  double s_coords[6] = {0, g->x, 0, g->y, 0, g->z};
  
  for(i=0; i<6; i++){
    if(p_angles[i] == 0){
      dist = D_INF;
    }
    else{
      dist = (s_coords[i] - p_coords[i])/p_angles[i];
      if(dist <= 0){
        dist = D_INF;
      }
    }
    if(dist < d){
      d = dist;
      g->surface_crossed = surfaces[i];
    }
  }

  return d;
}

// Returns the distance to the next collision for a particle
double distance_to_collision(Material *m)
{
  double d;

  if(m->xs_t == 0){
    d = D_INF;
  }
  else{
    d = -log(rn())/m->xs_t;
  }

  return d;
}

// Handles a particle crossing a surface in the geometry
void cross_surface(Particle *p, Geometry *g)
{
  // Handle vacuum boundary conditions (particle leaks out)
  if(g->bc == VACUUM){
    p->alive = FALSE;
  }

  // Handle reflective boundary conditions
  else if(g->bc == REFLECT){
    if(g->surface_crossed == X){
      p->u = -p->u;
      p->x = fabs(p->x);
    }
    else if(g->surface_crossed == Y){
      p->v = -p->v;
      p->y = fabs(p->y);
    }
    else if(g->surface_crossed == Z){
      p->w = -p->w;
      p->z = fabs(p->z);
    }
  }
  
  // Handle periodic boundary conditions
  else if(g->bc == PERIODIC){
    if(g->surface_crossed == X){
      p->x = g->x - p->x;
    }
    else if(g->surface_crossed == Y){
      p->y = g->y - p->y;
    }
    else if(g->surface_crossed == Z){
      p->z = g->z - p->z;
    }
  }

  return;
}

void collision(Particle *p, Material *m, Bank *fission_bank)
{
  int n;
  int i = 0;
  double prob = 0.0;
  double cutoff;
  Particle *p_new;
  Nuclide nuc;

  // Cutoff for sampling nuclide
  cutoff = rn()*m->xs_t;

  // Sample which nuclide particle has collision with
  while(prob < cutoff){
    nuc = m->nuclides[i];
    prob += nuc.atom_density * nuc.xs_t;
    i++;
  }

  // Cutoff for sampling reaction
  cutoff = rn()*nuc.xs_t;

  // Sample fission
  if(nuc.xs_f > cutoff){

    // Sample expected number of fission particles produced (arbitrary)
    n = rand() % 4;

    // Sample n new particles from the source distribution but at the current
    // particle's location
    for(i=0; i<n; i++){
      if(fission_bank->n >= fission_bank->sz){
        fission_bank->resize(fission_bank);
      }
      p_new = &(fission_bank->p[fission_bank->n]);
      p_new->alive = TRUE;
      p_new->energy = 1;
      p_new->last_energy = 0;
      p_new->mu = rn()*2 - 1;
      p_new->phi = rn()*2*PI;
      p_new->u = p_new->mu;
      p_new->v = sqrt(1 - p_new->mu*p_new->mu)*cos(p_new->phi);
      p_new->w = sqrt(1 - p_new->mu*p_new->mu)*sin(p_new->phi);
      p_new->x = p->x;
      p_new->y = p->y;
      p_new->z = p->z;
      fission_bank->n++;
    }
    p->alive = FALSE;
    p->event = FISSION;
  }

  // Sample absorption (disappearance)
  else if(nuc.xs_a > cutoff){
    p->alive = FALSE;
    p->event = ABSORPTION;
  }

  // Sample scattering
  else{
    p->mu = rn()*2 - 1;
    p->phi = rn()*2*PI;
    p->u = p->mu;
    p->v = sqrt(1 - p->mu*p->mu) * cos(p->phi);
    p->w = sqrt(1 - p->mu*p->mu) * sin(p->phi);
    p->event = SCATTER;
  }

  return;
}
