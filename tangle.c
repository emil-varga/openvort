#include <math.h>
#include "tangle.h"
#include "vortex_constants.h"

void alloc_arrays(struct tangle_state *tangle, size_t n)
{
  tangle->vnodes      = (struct vec3d*)malloc(sizeof(struct vec3d)*n);
  tangle->vnodes_new  = (struct vec3d*)malloc(sizeof(struct vec3d)*n);
  tangle->vels        = (struct vec3d*)malloc(sizeof(struct vec3d)*n);
  tangle->tangents    = (struct vec3d*)malloc(sizeof(struct vec3d)*n);
  tangle->normals     = (struct vec3d*)malloc(sizeof(struct vec3d)*n);

  tangle->connections = (struct neighbour_t*)malloc(sizeof(struct neighbour_t*)*n);

  tangle->N           = n;
  tangle->next_free   = 0;
  tangle->total_free  = n;
}

void free_arrays(struct tangle_state *tangle)
{
  free(tangle->vnodes);
  free(tangle->vnodes_new);
  free(tangle->vels);
  free(tangle->tangents);
  free(tangle->normals);
  free(tangle->connections);
}

void step_nodes(struct tangle_state *tangle, double dt)
{
  size_t i, k;
  for(i=0; i<tangle->N; ++i)
    for(k=0; k<3; ++k)
      tangle->vnodes[i].p[k] += dt*tangle->vels[i].p[k];
}

void update_tangent_normal(struct tangle_state *tangle, size_t k)
{
  struct vec3d *s0, *s1, *sm1;
  size_t i;

  //empty point has -1 connections
  if(tangle->connections[k].forward == -1 || tangle->connections[k].reverse == -1)
    return;

  s0  = &tangle->vnodes[k];
  s1  = &tangle->vnodes[tangle->connections[k].forward];
  sm1 = &tangle->vnodes[tangle->connections[k].reverse];

  double lf = vec3_dist(s0, s1);
  double lr = vec3_dist(s0, sm1);
  double d  = lf*lr*(lf + lr);

  for(i=0; i<3; ++i)
    {
      tangle->tangents[k].p[i] = s1->p[i]*lr*lr + s0->p[i]*(lf*lf + lr*lr) - sm1->p[i]*lf*lf;
      tangle->tangents[k].p[i] /= d;

      tangle->normals[k].p[i] = s1->p[i]*lr + sm1->p[i]*lf - s0->p[i]*(lf + lr);
      tangle->normals[k].p[i] /= d;
    }
}

struct vec3d segment_field(const struct tangle_state *tangle, size_t i, struct vec3d r)
{
  struct vec3d R;
  struct vec3d Rp1;

  vec3_sub(&R, tangle->vnodes+i, &r);
  vec3_sub(&Rp1, &tangle->vnodes[tangle->connections[i].forward], &r);

  double lR   = vec3_d(&R);
  double lRp1 = vec3_d(&Rp1);
  double denom = lR*lRp1*(lR*lRp1 + vec3_dot(&R, &Rp1));

  struct vec3d vv;

  vec3_cross(&vv, &R, &Rp1);
  vec3_mul(&vv, &vv, (lR + lRp1)/denom);

  return vv;
}

struct vec3d lia_velocity(const struct tangle_state *tangle, size_t i)
{
  const struct vec3d *p    = tangle->vnodes + i;
  const struct vec3d *next = tangle->vnodes + tangle->connections[i].forward;
  const struct vec3d *prev = tangle->vnodes + tangle->connections[i].reverse;

  double l_next = vec3_dist(p, next);
  double l_pref = vec3_dist(p, prev);

  double f = log(sqrt(l_next*l_pref)/VORTEX_WIDTH);

  struct vec3d vv;
  vec3_cross(&vv, tangle->tangents+i, tangle->normals+i);
  vec3_mul(&vv, &vv, f);

  return vv;
}

void update_velocity(struct tangle_state *tangle, size_t k)
{
  size_t m, i;

  tangle->vels[k] = lia_velocity(tangle, k);
  
  for(m=0; m<tangle->N; ++m)
    {
      if(tangle->connections[m].forward == -1 ||
	 m == k ||
	 m == tangle->connections[k].forward ||
	 m == tangle->connections[k].reverse)
	continue;
      
      struct vec3d segment_vel = segment_field(tangle, m, tangle->vnodes[k]);
      for(i=0; i<3; ++i)
	tangle->vels[k].p[i] += segment_vel.p[i];
    }
}

void update_velocities(struct tangle_state *tangle)
{
  size_t i;
  for(i=0; i<tangle->N; ++i)
    {
      update_velocity(tangle, i);
    }
}

void update_tangents_normals(struct tangle_state *tangle)
{
  size_t i;
  for(i=0; i<tangle->N; ++i)
    {
      update_tangent_normal(tangle, i);
    }
}

void update_tangle_next_free(struct tangle_state *tangle)
{
  tangle->next_free++;
}
