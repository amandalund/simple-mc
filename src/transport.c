#include "header.h"
#include "global.h"

// Main logic to move particle
void transport(Particle *p)
{
  double d_b;
  double d_c;
  double d;

  while(p->alive){

    // Recalculate macro xs if particle has changed energy
    if(p->energy != p->last_energy){
      calculate_xs();
    }

    // Find distance to boundary
    d_b = distance_to_boundary(p);

    // Find distance to collision
    d_c = distance_to_collision();

    // Take smaller of two distances
    d = d_b < d_c ? d_b : d_c;

    // Advance particle
    p->x = p->x + d*p->u;
    p->y = p->y + d*p->v;
    p->z = p->z + d*p->w;

    // Case where particle crosses boundary
    if(d_b < d_c){
      cross_surface(p);
    }
    // Case where particle has collision
    else{
      collision(p);

      // Score tallies
      if(tally->tallies_on == TRUE){
        score_tally(tally, p);
      }
    }
  }
  return;
}

// Calculates the macroscopic cross section of the material the particle is
// traveling through
void calculate_xs(void)
{
  int i;

  // Reset macroscopic cross sections to 0
  material->xs_t = 0.0;
  material->xs_f = 0.0;
  material->xs_a = 0.0;
  material->xs_s = 0.0;

  for(i=0; i<material->n_nuclides; i++){

    Nuclide nuc = material->nuclides[i];

    // Add contribution from this nuclide to total macro xs
    material->xs_t += nuc.atom_density * nuc.xs_t;

    // Add contribution from this nuclide to fission macro xs
    material->xs_f += nuc.atom_density * nuc.xs_f;

    // Add contribution from this nuclide to absorption macro xs
    material->xs_a += nuc.atom_density * nuc.xs_a;

    // Add contribution from this nuclide to scattering macro xs
    material->xs_s += nuc.atom_density * nuc.xs_s;
  }

  return;
}

// Returns the distance to the nearest boundary for a particle traveling in a
// certain direction
double distance_to_boundary(Particle *p)
{
  int i;
  double dist;
  double d = D_INF;
  int    surfaces[6] = {X0, X1, Y0, Y1, Z0, Z1};
  double p_angles[6] = {p->u, p->u, p->v, p->v, p->w, p->w};
  double p_coords[6] = {p->x, p->x, p->y, p->y, p->z, p->z};
  double s_coords[6] = {0, geometry->x, 0, geometry->y, 0, geometry->z};
  
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
      p->surface_crossed = surfaces[i];
    }
  }

  return d;
}

// Returns the distance to the next collision for a particle
double distance_to_collision()
{
  double d;

  if(material->xs_t == 0){
    d = D_INF;
  }
  else{
    d = -log(rn())/material->xs_t;
  }

  return d;
}

// Handles a particle crossing a surface in the geometry
void cross_surface(Particle *p)
{
  // Handle vacuum boundary conditions (particle leaks out)
  if(geometry->bc == VACUUM){
    p->alive = FALSE;
  }

  // Handle reflective boundary conditions
  else if(geometry->bc == REFLECT){
    if(p->surface_crossed == X0){
      p->u = -p->u;
      p->x = 0.0;
    }
    else if(p->surface_crossed == X1){
      p->u = -p->u;
      p->x = geometry->x;
    }
    else if(p->surface_crossed == Y0){
      p->v = -p->v;
      p->y = 0.0;
    }
    else if(p->surface_crossed == Y1){
      p->v = -p->v;
      p->y = geometry->y;
    }
    else if(p->surface_crossed == Z0){
      p->w = -p->w;
      p->z = 0.0;
    }
    else if(p->surface_crossed == Z1){
      p->w = -p->w;
      p->z = geometry->z;
    }
  }
  
  // Handle periodic boundary conditions
  else if(geometry->bc == PERIODIC){
    if(p->surface_crossed == X0){
      p->x = geometry->x;
    }
    else if(p->surface_crossed == X1){
      p->x = 0;
    }
    else if(p->surface_crossed == Y0){
      p->y = geometry->y;
    }
    else if(p->surface_crossed == Y1){
      p->y = 0;
    }
    else if(p->surface_crossed == Z0){
      p->z = geometry->z;
    }
    else if(p->surface_crossed == Z1){
      p->z = 0;
    }
  }

  return;
}

void collision(Particle *p)
{
  int n;
  int i = 0;
  double prob = 0.0;
  double cutoff;
  double nu = parameters->nu;
  Nuclide nuc = {0, 0, 0, 0, 0};

  // Cutoff for sampling nuclide
  cutoff = rn()*material->xs_t;

  // Sample which nuclide particle has collision with
  while(prob < cutoff){
    nuc = material->nuclides[i];
    prob += nuc.atom_density*nuc.xs_t;
    i++;
  }

  // Cutoff for sampling reaction
  cutoff = rn()*nuc.xs_t;

  // Sample fission
  if(nuc.xs_f > cutoff){

    // Sample number of fission neutrons produced
    if(rn() > nu - (int)nu){
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
      sample_fission_particle(&(fission_bank->p[fission_bank->n]), p);
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
