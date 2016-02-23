#include "simple_mc.h"

// Main logic to move particle
void transport(Parameters *parameters, Geometry *geometry, Material *material, Bank *source_bank, Bank *fission_bank, Tally *tally, int delay_tag, Queue *delay_queue, Particle *p)
{
  double d_b;
  double d_c;
  double d;

  while(p->alive){

    // Recalculate macro xs if particle has changed energy
    if(p->energy != p->last_energy){
      calculate_xs(material);
    }

    // Find distance to boundary
    d_b = distance_to_boundary(geometry, p);

    // Find distance to collision
    d_c = distance_to_collision(material);

    // Take smaller of two distances
    d = d_b < d_c ? d_b : d_c;

    // Advance particle
    p->x = p->x + d*p->u;
    p->y = p->y + d*p->v;
    p->z = p->z + d*p->w;

    // Case where particle crosses boundary
    if(d_b < d_c){
      cross_surface(geometry, p);
    }
    // Case where particle has collision
    else{
      collision(material, fission_bank, parameters->nu, delay_tag, delay_queue, p);

      // Score tallies
      if(tally->tallies_on == TRUE){
        score_tally(parameters, material, tally, p);
      }
    }
  }
  return;
}

// Calculates the macroscopic cross section of the material the particle is
// traveling through. Currently not used as the problem is one-group homogenous
// cube
void calculate_xs(Material *material)
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
double distance_to_boundary(Geometry *geometry, Particle *p)
{
  int i;
  double dist;
  double d = D_INF;
  int    surfaces[6] = {X0, X1, Y0, Y1, Z0, Z1};
  double p_angles[6] = {p->u, p->u, p->v, p->v, p->w, p->w};
  double p_coords[6] = {p->x, p->x, p->y, p->y, p->z, p->z};
  double s_coords[6] = {0, geometry->Lx, 0, geometry->Ly, 0, geometry->Lz};
  
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
double distance_to_collision(Material *material)
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
void cross_surface(Geometry *geometry, Particle *p)
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
      p->x = geometry->Lx;
    }
    else if(p->surface_crossed == Y0){
      p->v = -p->v;
      p->y = 0.0;
    }
    else if(p->surface_crossed == Y1){
      p->v = -p->v;
      p->y = geometry->Ly;
    }
    else if(p->surface_crossed == Z0){
      p->w = -p->w;
      p->z = 0.0;
    }
    else if(p->surface_crossed == Z1){
      p->w = -p->w;
      p->z = geometry->Lz;
    }
  }
  
  // Handle periodic boundary conditions
  else if(geometry->bc == PERIODIC){
    if(p->surface_crossed == X0){
      p->x = geometry->Lx;
    }
    else if(p->surface_crossed == X1){
      p->x = 0;
    }
    else if(p->surface_crossed == Y0){
      p->y = geometry->Ly;
    }
    else if(p->surface_crossed == Y1){
      p->y = 0;
    }
    else if(p->surface_crossed == Z0){
      p->z = geometry->Lz;
    }
    else if(p->surface_crossed == Z1){
      p->z = 0;
    }
  }

  return;
}

void collision(Material *material, Bank *fission_bank, double nu, int delay_tag, Queue *delay_queue, Particle *p)
{
  int nf;
  int i = 0;
  double prob = 0.0;
  double cutoff;
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
      nf = nu;
    }
    else{
      nf = nu + 1;
    }

    // Sample n new particles from the source distribution but at the current
    // particle's location. When building the delay bank, queue one particle to
    // delay bank and add all to fission bank. When using the delay bank, queue
    // one particle to delay bank, add remaining particles to fission bank, and
    // dequeue one particle from delay bank to add to fission bank
    if(fission_bank->n+nf >= fission_bank->sz){
      fission_bank->resize(fission_bank);
    }
    for(i=0; i<nf; i++){
      sample_fission_particle(&(fission_bank->p[fission_bank->n]), p);
      fission_bank->n++;
    }
    if(delay_tag != IGNORE_DB){
      enqueue(delay_queue, &(fission_bank->p[fission_bank->n-1]));
    }
    if(delay_tag == USE_DB){
      dequeue(delay_queue, &(fission_bank->p[fission_bank->n-1]));
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
