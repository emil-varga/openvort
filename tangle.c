#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "tangle.h"
#include "vortex_constants.h"
#include "util.h"

#define _DEBUG_

#ifdef _DEBUG_
#include <stdio.h>
#endif

#define _GNU_SOURCE
#include <fenv.h>

//this is for inserting new points
//fit an interpolating spline over 4 points and get a point
//somewhere in the middle later on
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

//struct to hold together the spline info
//one spline for each dimension
typedef struct _vspline {
  gsl_interp_accel *acc_xyz[3];
  gsl_spline *spline_xyz[3];
  double xi[4];
} vspline;

void free_vspline(vspline *spline)
{
  for(int k=0; k<3; ++k)
    {
      gsl_spline_free(spline->spline_xyz[k]);
      gsl_interp_accel_free(spline->acc_xyz[k]);
    }
  free(spline);
}

vspline *construct_local_spline(struct tangle_state *tangle,
				   int points[4])
{
  //3 coordinates for 4 points
  double xyz[3][4];

  //allocate splines
  vspline *spline = (vspline*)malloc(sizeof(vspline));
  for(int k=0; k<3; ++k) //loop over dimensions
    {
      spline->acc_xyz[k] = gsl_interp_accel_alloc();
      spline->spline_xyz[k] = gsl_spline_alloc(gsl_interp_cspline, 4);
    }
 
  //set up the data vectors
  for(int k=0; k<4; ++k) //loop over points
    {
      struct vec3d s = tangle->vnodes[points[k]];
      xyz[0][k] = s.p[0];
      xyz[1][k] = s.p[1];
      xyz[2][k] = s.p[2];
    }
  
  spline->xi[0] = 0;
  for(int k=1; k<4; ++k) //loop over points
    {
      spline->xi[k] = spline->xi[k-1] +
	vec3_dist(&tangle->vnodes[points[k]], &tangle->vnodes[points[k-1]]);
    }

  for(int k=0; k<3; ++k) //loop over dimensions
    gsl_spline_init(spline->spline_xyz[k], spline->xi, xyz[k], 4);

  return spline;
}

struct vec3d eval_vspline(vspline *spline, double where)
{
  struct vec3d out;

  for(int k=0; k<3; ++k) //loop over dimensions
    out.p[k] = gsl_spline_eval(spline->spline_xyz[k],
			       where, spline->acc_xyz[k]);

  return out;
}

//
// Implementation of public functions begins here
//
//

void alloc_arrays(struct tangle_state *tangle, size_t n)
{
  tangle->vnodes      = (struct vec3d*)malloc(sizeof(struct vec3d)*n);
  tangle->vnodes_new  = (struct vec3d*)malloc(sizeof(struct vec3d)*n);
  tangle->vels        = (struct vec3d*)malloc(sizeof(struct vec3d)*n);
  tangle->tangents    = (struct vec3d*)malloc(sizeof(struct vec3d)*n);
  tangle->normals     = (struct vec3d*)malloc(sizeof(struct vec3d)*n);

  tangle->connections = (struct neighbour_t*)malloc(n*sizeof(struct neighbour_t));
  for(int k=0; k<n; ++k)
    {
      tangle->connections[k].forward = -1;
      tangle->connections[k].reverse = -1;
    }

  tangle->N           = n;
  tangle->next_free   = 0;
  tangle->total_free  = n;
}

void expand_arrays(struct tangle_state *tangle, size_t n)
{
  int k;
  int old_n = tangle->N;
  tangle->vnodes      = (struct vec3d*)realloc(tangle->vnodes, sizeof(struct vec3d)*n);
  tangle->vnodes_new  = (struct vec3d*)realloc(tangle->vnodes_new, sizeof(struct vec3d)*n);
  tangle->vels        = (struct vec3d*)realloc(tangle->vels, sizeof(struct vec3d)*n);
  tangle->tangents    = (struct vec3d*)realloc(tangle->tangents, sizeof(struct vec3d)*n);
  tangle->normals     = (struct vec3d*)realloc(tangle->normals, sizeof(struct vec3d)*n);

  tangle->connections = (struct neighbour_t*)realloc(tangle->connections, n*sizeof(struct neighbour_t));

  if(old_n < n)
    {
      for(k=old_n; k<n; ++k)
  	{
  	  tangle->connections[k].forward = -1;
  	  tangle->connections[k].reverse = -1;
  	}
    }

  //only change the next free if the current next_free is not valid
  if(tangle->connections[tangle->next_free].forward == -1)
    tangle->next_free = tangle->N;
  tangle->N           = n;
  tangle->total_free  += n;
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
  };

  //first derivative
  //four point coefficients, nominators, O(d^4)
  double s_1_cf[] = {
    -d1*dm1*dm2,
    d2*dm1*dm2,
    -d1*d2*dm2,
    d1*d2*dm1
  };

  //second derivative
  //four point coefficients, nominators, O(d^3)
  double s_2_cf[] = {
    2*((dm1 - d1)*dm2 - d1*dm1),
    -2*((dm1 - d2)*dm2 - d2*dm1),
    2*((d2  + d1)*dm2 - d1*d2),
    -2*((d2  + d1)*dm1 - d1*d2)
  };
  
  for(i=0; i<3; ++i)
    {
      tangle->tangents[k].p[i] = 0;
      tangle->normals[k].p[i]  = 0;
      for(int z = 0; z<4; ++z)
	{
	  tangle->tangents[k].p[i] += s_1_cf[z]/d_s_diff[z]*ds[z].p[i];
	  tangle->normals[k].p[i]  += s_2_cf[z]/d_s_diff[z]*ds[z].p[i];
	}
    }
  double x = vec3_d(&tangle->tangents[k]);
  vec3_mul(&tangle->normals[k], &tangle->normals[k], 1/x/x);
  vec3_normalize(&tangle->tangents[k]);
}

static inline struct vec3d segment_field(const struct tangle_state *tangle, size_t i, struct vec3d r)
{
  struct vec3d R;
  struct vec3d Rp1;

  vec3_sub(&R, tangle->vnodes+i, &r);
  vec3_sub(&Rp1, &tangle->vnodes[tangle->connections[i].forward], &r);

  double lR   = vec3_d(&R);
  double lRp1 = vec3_d(&Rp1);
  double denom = lR*lRp1*(lR*lRp1 + vec3_dot(&R, &Rp1));
  double f = KAPPA/4/M_PI;

  struct vec3d vv;

  vec3_cross(&vv, &R, &Rp1);
  vec3_mul(&vv, &vv, f*(lR + lRp1)/denom);

  return vv;
}

static inline struct vec3d lia_velocity(const struct tangle_state *tangle, size_t i)
{
  const struct vec3d *p    = tangle->vnodes + i;
  const struct vec3d *next = tangle->vnodes + tangle->connections[i].forward;
  const struct vec3d *prev = tangle->vnodes + tangle->connections[i].reverse;

  double l_next = vec3_dist(p, next);
  double l_prev = vec3_dist(p, prev);

  double f = KAPPA*log(sqrt(l_next*l_prev)/VORTEX_WIDTH)/4/M_PI;

  struct vec3d vv;
  vec3_cross(&vv, tangle->tangents+i, tangle->normals+i);
  vec3_mul(&vv, &vv, f);

  return vv;
}

struct vec3d calculate_vs(struct tangle_state *tangle, struct vec3d r, size_t skip)
{
  int m;
  struct vec3d vs = vec3(0,0,0);

  for(m=0; m < tangle->N; ++m)
    {
      if(tangle->connections[m].forward == -1   ||
	 skip == m                              ||
	 skip == tangle->connections[m].forward ||
	 skip == tangle->connections[m].reverse)
	continue;
      struct vec3d ivs = segment_field(tangle, m, r);
      vec3_add(&vs, &vs, &ivs);
    }

  return vs;
}

void update_velocity(struct tangle_state *tangle, size_t k)
{
  size_t m, i;
  if(tangle->connections[k].forward == -1)
    return;

  tangle->vels[k] = lia_velocity(tangle, k);
  //tangle->vels[k] = vec3(0,0,0);
  
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

static inline int search_next_free(struct tangle_state *tangle)
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

int get_tangle_next_free(struct tangle_state *tangle)
{
  int idx = search_next_free(tangle);

  if(idx < 0)
    {
      expand_arrays(tangle, 2*tangle->N);
      idx = search_next_free(tangle);
    }

  return idx;
}

int num_free_points(struct tangle_state *tangle)
{
  int sum=0;
  for(int k=0; k<tangle->N; ++k)
    if(tangle->connections[k].forward < 0)
      sum++;
  return sum;
}

void remove_point(struct tangle_state *tangle, int point_idx);
void add_point(struct tangle_state *tangle, int point_idx);
void remesh(struct tangle_state *tangle, double min_dist, double max_dist)
{
  for(int k=0; k<tangle->N; ++k)
    {
      if(tangle->connections[k].forward < 0) //empty point
	continue;

      int next = tangle->connections[k].forward;
      int prev = tangle->connections[k].reverse;
  
      double lf = vec3_dist(tangle->vnodes + k,
			    tangle->vnodes + next);
      double lr = vec3_dist(tangle->vnodes + k,
			    tangle->vnodes + prev);

      //can we remove point k?
      if( (lf < min_dist || lr < min_dist) && (lf + lr) < max_dist )
	  remove_point(tangle, k);

      //do we need an extra point?
      if( lf > max_dist ) //since we are adding between k and next, check only lf
	add_point(tangle, k);
    }
}

void remove_point(struct tangle_state *tangle, int point_idx)
{
  int prev = tangle->connections[point_idx].reverse;
  int next = tangle->connections[point_idx].forward;
  tangle->connections[prev].forward = next;
  tangle->connections[next].reverse = prev;
  tangle->connections[point_idx].reverse =
    tangle->connections[point_idx].forward = -1;
}


//add a point between p and p+1 (p+1 in the sense of connections)
void add_point(struct tangle_state *tangle, int p)
{
  int next = tangle->connections[p].forward;
  int new_pt = get_tangle_next_free(tangle);

  struct vec3d s0 = tangle->vnodes[p];
  struct vec3d s1 = tangle->vnodes[next];

  struct vec3d s0p = tangle->tangents[p];
  struct vec3d s1p = tangle->tangents[next];

  struct vec3d s0pp = tangle->normals[p];
  struct vec3d s1pp = tangle->normals[next];

  double l = vec3_dist(&s0, &s1);

  struct vec3d a, b, c, new;
  vec3_add(&a, &s0, &s1);
  vec3_mul(&a, &a, 0.5);

  vec3_sub(&b, &s0p, &s1p);
  vec3_mul(&b, &b, l/4);

  vec3_add(&c, &s0pp, &s1pp);
  vec3_mul(&c, &c, l*l/16);

  vec3_add(&new, &a, &b);
  vec3_add(&new, &new, &c);

  tangle->vnodes[new_pt] = new;
  tangle->connections[new_pt].reverse = p;
  tangle->connections[new_pt].forward = next;
  tangle->connections[p].forward = new_pt;
  tangle->connections[next].reverse = new_pt;
}
