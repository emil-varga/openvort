#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>
#include <assert.h>
#include "tangle.h"
#include "vortex_constants.h"
#include "normal_fluid.h"
#include "util.h"
#include "vec3_maths.h"

#ifdef _DEBUG_
#include <stdio.h>
#endif

#define _GNU_SOURCE
#include <fenv.h>


//
// Implementation of public functions begins here
//
//

/*
 * Inward-facing normals of the box boundary face walls
 */
const struct vec3d boundary_normals[] = {
    {{1, 0, 0}},
    {{-1, 0, 0}},
    {{0, 1, 0}},
    {{0, -1, 0}},
    {{0, 0, 1}},
    {{0, 0, -1}}
};


void create_tangle(struct tangle_state *tangle, size_t n)
{
  tangle->vnodes      = (struct vec3d*)malloc(sizeof(struct vec3d)*n);
  tangle->vnodes_new  = (struct vec3d*)malloc(sizeof(struct vec3d)*n);
  tangle->vels        = (struct vec3d*)malloc(sizeof(struct vec3d)*n);
  tangle->vs          = (struct vec3d*)malloc(sizeof(struct vec3d)*n);
  tangle->tangents    = (struct vec3d*)malloc(sizeof(struct vec3d)*n);
  tangle->normals     = (struct vec3d*)malloc(sizeof(struct vec3d)*n);
  tangle->recalculate = (int*)malloc(sizeof(int)*n);
  tangle->status      = (node_status*)malloc(sizeof(node_status)*n);

  tangle->connections = (struct neighbour_t*)malloc(n*sizeof(struct neighbour_t));
  for(size_t k=0; k<n; ++k)
    {
      tangle->status[k].status = EMPTY;
      tangle->connections[k].forward = -1;
      tangle->connections[k].reverse = -1;
    }

  //default initialisation is to open bounadry conditions
  struct domain_box box = {
      .bottom_left_front = {{0,0,0}},
      .top_right_back = {{1,1,1}},
      .wall = {WALL_OPEN, WALL_OPEN, WALL_OPEN, WALL_OPEN, WALL_OPEN, WALL_OPEN}
  };

  tangle->box = box;

  tangle->N           = n;
  tangle->next_free   = 0;
  tangle->total_free  = n;
}

void expand_tangle(struct tangle_state *tangle, size_t n)
{
  size_t k;
  size_t old_n = tangle->N;
  tangle->vnodes      = (struct vec3d*)realloc(tangle->vnodes, sizeof(struct vec3d)*n);
  tangle->vnodes_new  = (struct vec3d*)realloc(tangle->vnodes_new, sizeof(struct vec3d)*n);
  tangle->vels        = (struct vec3d*)realloc(tangle->vels, sizeof(struct vec3d)*n);
  tangle->vs          = (struct vec3d*)realloc(tangle->vs, sizeof(struct vec3d)*n);
  tangle->tangents    = (struct vec3d*)realloc(tangle->tangents, sizeof(struct vec3d)*n);
  tangle->normals     = (struct vec3d*)realloc(tangle->normals, sizeof(struct vec3d)*n);
  tangle->recalculate = (int*)realloc(tangle->recalculate, sizeof(int)*n);
  tangle->status      = (node_status*)realloc(tangle->status, sizeof(node_status)*n);

  tangle->connections = (struct neighbour_t*)realloc(tangle->connections, n*sizeof(struct neighbour_t));

  if(old_n < n)
    {
      for(k=old_n; k<n; ++k)
  	{
	  tangle->status[k].status = EMPTY;
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

void free_tangle(struct tangle_state *tangle)
{
  free(tangle->vnodes);
  free(tangle->vnodes_new);
  free(tangle->vels);
  free(tangle->vs);
  free(tangle->tangents);
  free(tangle->normals);
  free(tangle->connections);
  free(tangle->recalculate);
  free(tangle->status);
}

struct vec3d step_node(struct tangle_state *tangle, int i, int where)
{
  assert(tangle->status[i].status != EMPTY);

  if(where == 0)
    return tangle->vnodes[i];

  if(tangle->status[i].status == FREE)
    {
      if(where > 0)
	return step_node(tangle, tangle->connections[i].forward, where-1);
      else if (where < 0)
	return step_node(tangle, tangle->connections[i].reverse, where+1);
    }
  else if(tangle->status[i].status == PINNED ||
          tangle->status[i].status == PINNED_SLIP)
    {
      struct vec3d out;
      if(where > 0)
	{
	  if(tangle->connections[i].forward < 0)
	    out = step_node(tangle, i, -where);
	  else
	    return step_node(tangle, tangle->connections[i].forward, where-1);
	}
      else if(where < 0)
	{
	  if(tangle->connections[i].reverse < 0)
	    out = step_node(tangle, i, -where);
	  else
	    return step_node(tangle, tangle->connections[i].reverse, where+1);
	}

      //if we are here it means we have ran into a wall and we need to flip the node
      return mirror_shift(&out, &tangle->box, tangle->status[i].pin_wall);
    }
  else
    {
      error("Walking across empty node.");//we should never get here
      return vec3(0,0,0);
    }
}

void update_tangent_normal(struct tangle_state *tangle, size_t k)
{
  struct vec3d s0, s1, sm1;
  struct vec3d s2, sm2;
  size_t i;

  //vector differences
  struct vec3d ds[4];
  struct segment dseg[4];
  struct segment dseg_12;
  struct segment dseg_m12;

  //empty point has -1 connections
  if(tangle->status[k].status == EMPTY)
    return;

  int next = tangle->connections[k].forward;
  int prev = tangle->connections[k].reverse;
  int next2, prev2;


  next2 = tangle->connections[next].forward;
  prev2 = tangle->connections[prev].reverse;

  /*
   * TODO: sm2..s2 should be constructed iteratively based on pinning and
   * periodic conditions using shift and mirror functions from vec3_maths
   */
  s0  = tangle->vnodes[k];
  s1 = step_node(tangle, k, 1);
  s2 = step_node(tangle, k, 2);
  sm1 = step_node(tangle, k, -1);
  sm2 = step_node(tangle, k, -2);

  dseg[0] = seg_pwrap(&s0, &s2, &tangle->box);
  dseg[1] = seg_pwrap(&s0, &s1, &tangle->box);
  dseg[2] = seg_pwrap(&s0, &sm1, &tangle->box);
  dseg[3] = seg_pwrap(&s0, &sm2, &tangle->box);
  dseg_12 = seg_pwrap(&s1, &s2, &tangle->box);
  dseg_m12 = seg_pwrap(&sm1, &sm2, &tangle->box);

  for(int j = 0; j<4; ++j)
    ds[j] = segment_to_vec(&dseg[j]);

  double d1 = segment_len(&dseg[1]);
  double d2 = d1 + segment_len(&dseg_12);
  double dm1 = segment_len(&dseg[2]);
  double dm2 = dm1 + segment_len(&dseg_m12);

//  double d1 = vec3_dist(&s0, &s1);
//  double d2 = d1 + vec3_dist(&s1, &s2);
//
//  double dm1 = vec3_dist(&s0, &sm1);
//  double dm2 = dm1 + vec3_dist(&sm1, &sm2);

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

/*
 * Calculates the velocity field of a straight vortex segment given by seg at position r
 */
static inline struct vec3d segment_field1(struct segment *seg, struct vec3d r)
{
  struct vec3d R;
  struct vec3d Rp1;

  vec3_sub(&R, &seg->r1, &r);
  vec3_sub(&Rp1, &seg->r2, &r);

  double lR   = vec3_d(&R);
  double lRp1 = vec3_d(&Rp1);
  double denom = lR*lRp1*(lR*lRp1 + vec3_dot(&R, &Rp1));
  double f = KAPPA/4/M_PI;

  struct vec3d vv;

  vec3_cross(&vv, &R, &Rp1);
  vec3_mul(&vv, &vv, f*(lR + lRp1)/denom);

  return vv;
}

/*
 * Calculates the velocity field of a segment in a given tangle
 */
static inline struct vec3d segment_field(const struct tangle_state *tangle, size_t i, struct vec3d r)
{
  struct segment seg = {
      .r1 = tangle->vnodes[i],
      .r2 = tangle->vnodes[tangle->connections[i].forward]
  };
  //TODO: this should create the segment in accordance with boundary conditions

  return segment_field1(&seg, r);
}

static inline struct vec3d lia_velocity(const struct tangle_state *tangle, size_t i)
{
  //TODO: this needs to be aware of box boundary conditions
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

struct vec3d calculate_vs_shift(struct tangle_state *tangle, struct vec3d r, int skip,
				const struct vec3d *shift)
{
  int m;
  struct vec3d vs = vec3(0,0,0);

  if(shift)
    vec3_add(&r, &r, shift);


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

struct vec3d calculate_vs(struct tangle_state *tangle, struct vec3d r, int skip)
{
  return calculate_vs_shift(tangle, r, skip, NULL);
}

void update_velocity(struct tangle_state *tangle, int k)
{
  int m, i;
  if(tangle->connections[k].forward == -1)
    return;

  if(tangle->status[k].status == PINNED)
    {
      tangle->vs[k] = vec3(0,0,0);
      return;
    }

  tangle->vs[k] = lia_velocity(tangle, k);
  
  for(m=0; m<tangle->N; ++m)
    {
      if(tangle->connections[m].forward == -1 ||
  	 m == k                               ||
  	 m == tangle->connections[k].forward  ||
  	 m == tangle->connections[k].reverse)
  	continue;
      
      struct vec3d segment_vel = segment_field(tangle, m, tangle->vnodes[k]);
      for(i=0; i<3; ++i)
  	tangle->vs[k].p[i] += segment_vel.p[i];
    }

  /*
   * TODO: here, calculate_vs_shifts should be called with appropriate shifts
   * for the boundary conditions and also the appropriate sign changes for mirrors
   */

  tangle->vels[k] = tangle->vs[k];

  if(use_mutual_friction)
    {
      struct vec3d tmp, dv;

      //the velocity difference
      get_vn(&tangle->vnodes[k], &dv);
      vec3_sub(&dv, &dv, &tangle->vs[k]);

      //the dissipative term
      vec3_cross(&tmp, &tangle->tangents[k], &dv);
      vec3_mul(&tmp, &tmp, alpha);
      vec3_add(&tangle->vels[k], &tangle->vels[k], &tmp);

      //the non-dissipative term
      vec3_cross(&tmp, &tangle->tangents[k], &dv);
      vec3_cross(&tmp, &tangle->tangents[k], &tmp);
      vec3_mul(&tmp, &tmp, -alpha_p);
      vec3_add(&tangle->vels[k], &tangle->vels[k], &tmp);
    }

  if(tangle->status[k].status == PINNED_SLIP)
    {
      //remove the component of velocity normal to the wall
    }
}

void update_velocities(struct tangle_state *tangle)
{
  int i;
  #pragma omp parallel private(i) num_threads(6)
    {
      #pragma omp for
      for(i=0; i<tangle->N; ++i)
	{
	  update_velocity(tangle, i);
	}
    }
}

void update_tangents_normals(struct tangle_state *tangle)
{
  int i;
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
      expand_tangle(tangle, 2*tangle->N);
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

static inline int out_of_box(const struct vec3d *where)
{
  double x = where->p[0];
  double y = where->p[1];
  double z = where->p[2];
  if(x < computation_box[X_L])
    return X_L;
  if(x > computation_box[X_H])
    return X_H;
  if(y < computation_box[Y_L])
    return Y_L;
  if(y > computation_box[Y_H])
    return Y_H;
  if(z < computation_box[Z_L])
    return Z_L;
  if(z > computation_box[Z_H])
    return Z_H;

  return 0;
}
void enforce_boundaries(struct tangle_state *tangle)
{
  int face;
  for(int k=0; k < tangle->N; ++k)
    {
      if(tangle->status[k].status == EMPTY)
	continue;
      face = out_of_box(&tangle->vnodes[k]);
      if(face)
	  {
	    assert_msg(tangle->status[k].status != PINNED ||
		       tangle->status[k].status != PINNED_SLIP,
		       "pinned node outside of the box\n"
		       "this should have been caut with reconnections")
	    switch(face)
	    {
	      case X_L:
		break;
	      default:
		assert_msg(0, "This is not a valid face.");
	    }
	  }
    }
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

void eliminate_small_loops(struct tangle_state *tangle, int loop_length)
{
  for(int k=0; k < tangle->N; ++k)
    tangle->recalculate[k] = 0;

  for(int k=0; k < tangle->N; ++k)
    {
      if(tangle->connections[k].forward < 0 ||
	 tangle->recalculate[k])
	continue; //empty or visited point

      tangle->recalculate[k]++;

      int loop = 0;
      int z = tangle->connections[k].forward;
      while(z != k)
	{
	  tangle->recalculate[z]++;
	  z = tangle->connections[z].forward;
	  loop++;
	}
      if(loop < loop_length)
	{ //the loop is short, delete it
	  #ifdef _DEBUG_
	  printf("Eliminating loop starting at %d\n", k);
	  #endif
	  z = k;
	  while(tangle->connections[z].forward > 0)
	    {
	      int tmp = z;
	      z = tangle->connections[z].forward;
	      tangle->connections[tmp].forward = -1;
	      tangle->connections[tmp].reverse = -1;
	      tangle->status[tmp].status = EMPTY;
	    }
	}
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

  tangle->status[point_idx].status = EMPTY;
}


//add a point between p and p+1 (p+1 in the sense of connections)
void add_point(struct tangle_state *tangle, int p)
{
  int next = tangle->connections[p].forward;
  int new_pt = get_tangle_next_free(tangle);
  tangle->status[new_pt].status = FREE;

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
