#include <math.h>
#include <string.h>
#include "tangle.h"
#include "vortex_constants.h"

void alloc_arrays(struct tangle_state *tangle, size_t n)
{
  tangle->vnodes      = (struct vec3d*)malloc(sizeof(struct vec3d)*n);
  tangle->vnodes_new  = (struct vec3d*)malloc(sizeof(struct vec3d)*n);
  tangle->vels        = (struct vec3d*)malloc(sizeof(struct vec3d)*n);
  tangle->tangents    = (struct vec3d*)malloc(sizeof(struct vec3d)*n);
  tangle->normals     = (struct vec3d*)malloc(sizeof(struct vec3d)*n);

  tangle->connections = (struct neighbour_t*)malloc(n*sizeof(struct neighbour_t*));
  for(int k=0; k<n; ++k)
    {
      tangle->connections[k].forward = -1;
      tangle->connections[k].reverse = -1;
    }

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

void step_nodes2(struct tangle_state *result,
		 const struct tangle_state *tangle, double dt)
{
  size_t i, k;
  if(tangle != result)
    {
      //copy necessary things from tangle to result
      memcpy(result->connections, tangle->connections, tangle->N);

      //velocities and tangents/normals will be recalculated after stepping
    }
  for(i=0; i<tangle->N; ++i)
    for(k=0; k<3; ++k)
      result->vnodes[i].p[k] += dt*tangle->vels[i].p[k];

  //calculate the normals/tangents for the new tangle
  update_tangle(result);
}

void update_tangent_normal(struct tangle_state *tangle, size_t k)
{
  struct vec3d *s0, *s1, *sm1;
  struct vec3d *s2, *sm2;
  size_t i;

  //vector differences
  struct vec3d ds[4];

  //empty point has -1 connections
  if(tangle->connections[k].forward == -1 ||
     tangle->connections[k].reverse == -1)
    return;

  int next = tangle->connections[k].forward;
  int prev = tangle->connections[k].reverse;

  s0  = tangle->vnodes + k;
  s1  = tangle->vnodes + next;
  sm1 = tangle->vnodes + prev;

  s2 = tangle->vnodes + tangle->connections[next].forward;
  sm2 = tangle->vnodes + tangle->connections[prev].reverse;

  vec3_sub(ds + 0, s2, s0);
  vec3_sub(ds + 1, s1, s0);
  vec3_sub(ds + 2, sm1, s0);
  vec3_sub(ds + 3, sm2, s0);

  double d1 = vec3_dist(s0, s1);
  double d2 = d1 + vec3_dist(s1, s2);

  double dm1 = vec3_dist(s0, sm1);
  double dm2 = dm1 + vec3_dist(sm1, sm2);

  //four point coefficients, denominators
  double d_s_diff[] = {
    d2*(d2 - d1)*(dm1 + d2)*(dm2 + d2),
    d1*(d2 - d1)*(dm1 + d1)*(dm2 + d1),
    dm1*(dm1 + d1)*(dm1 + d2)*(dm2 - dm1),
    dm2*(dm2 + d1)*(dm2 + d2)*(dm2 - dm1)
  }

  //first derivative
  //four point coefficients, nominators, O(d^4)
  double s_1_cf[] = {
    -d1*dm1*dm2,
    d2*dm1*dm2,
    -d1*d2*dm2,
    d1*d2*dm1
  }

  //second derivative
  //four point coefficients, nominators, O(d^3)
  double s_2_cf[] = {
     2*((dm1 - d1)*dm2 - d1*dm1),
    -2*((dm1 - d2)*dm2 - d2*dm1),
     2*((d2  + d1)*dm2 - d1*d2),
    -2*((d2  + d1)*dm1 - d1*d2)
  }
  
  for(i=0; i<3; ++i)
    {
      for(int z = 0; z<4; ++z)
	{
	  tangle->tangents[k].p[i] = s_1_cf[z]/d_s_diff[z]*ds[z][i];
	  tangle->normals[k].p[i]  = s_2_cf[z]/d_s_diff[z]*ds[z][i];
	}
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
  double l_prev = vec3_dist(p, prev);

  double f = log(sqrt(l_next*l_prev)/VORTEX_WIDTH);

  struct vec3d vv;
  vec3_cross(&vv, tangle->tangents+i, tangle->normals+i);
  vec3_mul(&vv, &vv, f);

  return vv;
}

void update_velocity(struct tangle_state *tangle, size_t k)
{
  size_t m, i;
  if(tangle->connections[k].forward == -1)
    return;

  tangle->vels[k] = lia_velocity(tangle, k);
  
  for(m=0; m<tangle->N; ++m)
    {
      if(tangle->connections[m].forward == -1 ||
	 m == k                               ||
	 m == tangle->connections[k].forward  ||
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

int get_tangle_next_free(struct tangle_state *tangle)
{
  //if something is already available, just return that
  if(tangle->next_free < tangle->N &&
     tangle->connections[tangle->next_free].forward == -1)
    return tangle->next_free++;

  //otherwise we have to search for it

  for(int k=0; k<tangle->N; ++k)
    if(tangle->connections[k].forward == -1)
      {
	tangle->next_free = k+1;
	return k;
      }

  //we haven't found anything, return error
  return -1;
}

int num_free_points(struct tangle_state *tangle)
{
  int sum=0;
  for(int k=0; k<tangle->N; ++k)
    if(tangle->connections[k].forward < 0)
      sum++;
  return sum;
}

void remesh(struct tangle_state *tangle,
	    double min_dist, double max_dist)
{
  for(int k=0; k<tangle->N; ++k)
    {
      if(tangle->connections[k].forward < 0) //empty point
	continue;

      int next = tangle->connections[k].forward;
      int prev = tangle->connections[k].reverse;

      double lf = vec3_dist(tangle->vnodes + k
			    tangle->vnodes + next);
      double lr = vec3_dist(tangle->vnodes + k
			    tangle->vnodes + prev);

      //can we remove point k?
      if( (lf < min_dist || lr < min_dist) && (lf + lr) < max_dist )
	{
	  tangle->connections[prev].forward = next;
	  tangle->connections[next].reverse = prev;
	  tangle->connections[k].forward = -1;
	  tangle->connections[k].reverse = -1;
	}

      //do we need an extra point?
      if( lf > max_dist ) //since we are adding between k and next, check only lf
	{
	  int new_pt = get_tangle_next_free(tangle);
	}
    }
}
