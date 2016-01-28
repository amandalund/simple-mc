#include "header.h"

// Main logic to move particle
void transport(Particle *p, Geometry *g, Material *m, Tally *t, Bank *fission_bank, double keff, Parameters *params)
{
  while(p->alive){

    // Recalculate macro xs if particle has changed energy
    if(p->energy != p->last_energy){
      calculate_xs(p, m);
    }

    // Find distance to boundary
    double d_b = distance_to_boundary(p, g);

    // Find distance to collision
    double d_c = distance_to_collision(m, params);

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
      collision(p, m, fission_bank, keff, params);

      // Score tallies
      if(t->tallies_on == TRUE){
        score_tally(t, p, m, params);
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
  m->xs_s = 0.0;

  for(i=0; i<m->n_nuclides; i++){

    Nuclide nuc = m->nuclides[i];

    // Add contribution from this nuclide to total macro xs
    m->xs_t += nuc.atom_density * nuc.xs_t;

    // Add contribution from this nuclide to fission macro xs
    m->xs_f += nuc.atom_density * nuc.xs_f;

    // Add contribution from this nuclide to absorption macro xs
    m->xs_a += nuc.atom_density * nuc.xs_a;

    // Add contribution from this nuclide to scattering macro xs
    m->xs_s += nuc.atom_density * nuc.xs_s;
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
  int    surfaces[6] = {X0, X1, Y0, Y1, Z0, Z1};
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
double distance_to_collision(Material *m, Parameters *params)
{
  double d;

  if(m->xs_t == 0){
    d = D_INF;
  }
  else{
    d = -log(rn(&(params->seed)))/m->xs_t;
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
    if(g->surface_crossed == X0){
      p->u = -p->u;
      p->x = 0.0;
    }
    else if(g->surface_crossed == X1){
      p->u = -p->u;
      p->x = g->x;
    }
    else if(g->surface_crossed == Y0){
      p->v = -p->v;
      p->y = 0.0;
    }
    else if(g->surface_crossed == Y1){
      p->v = -p->v;
      p->y = g->y;
    }
    else if(g->surface_crossed == Z0){
      p->w = -p->w;
      p->z = 0.0;
    }
    else if(g->surface_crossed == Z1){
      p->w = -p->w;
      p->z = g->z;
    }
  }
  
  // Handle periodic boundary conditions
  else if(g->bc == PERIODIC){
    if(g->surface_crossed == X0){
      p->x = g->x;
    }
    else if(g->surface_crossed == X1){
      p->x = 0;
    }
    else if(g->surface_crossed == Y0){
      p->y = g->y;
    }
    else if(g->surface_crossed == Y1){
      p->y = 0;
    }
    else if(g->surface_crossed == Z0){
      p->z = g->z;
    }
    else if(g->surface_crossed == Z1){
      p->z = 0;
    }
  }

  return;
}

void collision(Particle *p, Material *m, Bank *fission_bank, double keff, Parameters *params)
{
  int n;
  int i = 0;
  double prob = 0.0;
  double cutoff;
  double nu = params->nu;
  Nuclide nuc;

  // Cutoff for sampling nuclide
  cutoff = rn(&(params->seed))*m->xs_t;

  // Sample which nuclide particle has collision with
  while(prob < cutoff){
    nuc = m->nuclides[i];
    prob += nuc.atom_density*nuc.xs_t;
    i++;
  }

  // Cutoff for sampling reaction
  cutoff = rn(&(params->seed))*nuc.xs_t;

  // Sample fission
  if(nuc.xs_f > cutoff){

    // Sample number of fission neutrons produced
    if(rn(&(params->seed)) > nu - (int)nu){
      n = nu;
    }
    else{
      n = nu + 1;
    }

    // Sample n new particles from the source distribution but at the current
    // particle's location
    if(fission_bank->n+n >= fission_bank->sz){
      fission_bank->resize(fission_bank);
    }
    for(i=0; i<n; i++){
      sample_fission_particle(&(fission_bank->p[fission_bank->n]), p, params);
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
    p->mu = rn(&(params->seed))*2 - 1;
    p->phi = rn(&(params->seed))*2*PI;
    p->u = p->mu;
    p->v = sqrt(1 - p->mu*p->mu) * cos(p->phi);
    p->w = sqrt(1 - p->mu*p->mu) * sin(p->phi);
    p->event = SCATTER;
  }

  return;
}
