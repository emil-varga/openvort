/*******************************************************************************
 * Copyright (C) 2018 Emil Varga <varga.emil@gmail.com>
 * 
 * This file is part of OpenVort
 * 
 * OpenVort is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>
#include <assert.h>
#include <external_velocity.h>
#include "tangle.h"
#include "vortex_constants.h"
#include "util.h"
#include "vec3_maths.h"
#include "octree.h"

#ifdef _DEBUG_
#include <stdio.h>
#endif

#define _GNU_SOURCE
#include <fenv.h>

struct vec3d shifted(const struct image_tangle *shift, const struct tangle_state *tangle,
		     const struct vec3d *r)
{
  /*
   * Mirror images of only depth 1 are supported (higher depth could make sense
   * for two mirror facing each other)
   */
  struct vec3d rs = *r;
  double Ls[3];
  for(int k = 0; k<3; ++k)
    Ls[k] = tangle->box.top_right_front.p[k] - tangle->box.bottom_left_back.p[k];

  for(int k=0; k<3; ++k)
    {
      rs.p[k] -= Ls[k]*shift->shift[k];
    }

  switch(shift->reflect)
  {
    case X_L:
      rs.p[0] = 2*tangle->box.bottom_left_back.p[0] - r->p[0];
      break;
    case X_H:
      rs.p[0] = 2*tangle->box.top_right_front.p[0] - r->p[0];
      break;
    case Y_L:
      rs.p[1] = 2*tangle->box.bottom_left_back.p[1] - r->p[1];
      break;
    case Y_H:
      rs.p[1] = 2*tangle->box.top_right_front.p[1] - r->p[1];
      break;
    case Z_L:
      rs.p[2] = 2*tangle->box.bottom_left_back.p[2] - r->p[2];
      break;
    case Z_H:
      rs.p[2] = 2*tangle->box.top_right_front.p[2] - r->p[2];
      break;
    default: //no reflection
      break;
  }

  return rs;
}


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
      tangle->status[k].pin_wall = -1;
      tangle->connections[k].forward = -1;
      tangle->connections[k].reverse = -1;
    }

  //default initialisation is to open bounadry conditions
  struct domain_box box = {
      .bottom_left_back = {{0,0,0}},
      .top_right_front = {{1,1,1}},
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
	  tangle->status[k].pin_wall = -1;
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

struct vec3d step_node(const struct tangle_state *tangle, int i, int where)
{
  assert(tangle->status[i].status != EMPTY);

  if(where == 0)
    return tangle->vnodes[i];

  if(tangle->status[i].status == FREE) {
    if(where > 0)
	    return step_node(tangle, tangle->connections[i].forward, where-1);
    else if (where < 0)
	    return step_node(tangle, tangle->connections[i].reverse, where+1);
  } else if(tangle->status[i].status == PINNED || tangle->status[i].status == PINNED_SLIP) {
    struct vec3d out;
    if(where > 0) {
	    if(tangle->connections[i].forward < 0)
	      out = step_node(tangle, i, -where);
	    else
	      return step_node(tangle, tangle->connections[i].forward, where-1);
	  } else if(where < 0) {
	    if(tangle->connections[i].reverse < 0)
	      out = step_node(tangle, i, -where);
	    else
	      return step_node(tangle, tangle->connections[i].reverse, where+1);
	  }
      //if we are here it means we have ran into a wall and we need to flip the node
      return mirror_shift(&out, &tangle->box, tangle->status[i].pin_wall);
    } else {
      error("Walking across empty node.");//we should never get here
      return vec3(0,0,0);
    }

  //to suppress warning
  return vec3(0,0,0);
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

  if(tangle->status[k].status == EMPTY)
    return;

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
  //double x = vec3_d(&tangle->tangents[k]);
  //vec3_mul(&tangle->normals[k], &tangle->normals[k], 1/x/x);
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

  //this can happen in periodic boundary conditions
  //TODO: the logic should be moved higher
  if(lR < 1e-8 || lRp1 < 1e-8)
    return vec3(0, 0, 0);

  //if R and Rp1 are colinear, the result is 0
  //but code below would try to calculate 0/0
  if(fabs(vec3_ndot(&R, &Rp1) - 1) < 1e-8)
    return vec3(0,0,0);

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
  int next = tangle->connections[i].forward;
  if(next == -1) //this is an edge point on a wall
    return vec3(0,0,0);
  struct segment seg = seg_pwrap(tangle->vnodes + i, tangle->vnodes + next, &tangle->box);

  return segment_field1(&seg, r);
}

static inline struct vec3d lia_velocity(const struct tangle_state *tangle, int i)
{
  const struct vec3d *p    = tangle->vnodes + i;
  struct vec3d next = step_node(tangle, i, +1);
  struct vec3d prev = step_node(tangle, i, -1);

  struct segment sf = seg_pwrap(p, &next, &tangle->box);
  struct segment sr = seg_pwrap(&prev, p, &tangle->box);

  double l_next = segment_len(&sf);
  double l_prev = segment_len(&sr);

  double f = KAPPA*log(2*sqrt(l_next*l_prev)/sqrt(sqrt(M_E))/VORTEX_WIDTH)/4/M_PI;

  struct vec3d vv;
  vec3_cross(&vv, &tangle->tangents[i], &tangle->normals[i]);
  vec3_mul(&vv, &vv, f);

  return vv;
}

struct vec3d calculate_vs_shift(struct tangle_state *tangle, struct vec3d r, int skip, const struct vec3d *shift)
{
  int m;
  struct vec3d vs = vec3(0,0,0);

  if(shift)
    vec3_add(&r, &r, shift);


  for(m=0; m < tangle->N; ++m) {
    if(tangle->connections[m].forward == -1   ||
	     skip == m                              ||
	     skip == tangle->connections[m].forward)
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

void update_velocity(struct tangle_state *tangle, int k, double t)
{
  int m, i;
  if(tangle->status[k].status == EMPTY)
    return;

  if(tangle->status[k].status == PINNED)
    {
      tangle->vs[k] = vec3(0,0,0);
      //get the boundary velocity, by default non-moving
      get_vb(&tangle->vnodes[k], t, &tangle->vels[k]);
      return;
    }

  //save the LIA velocity for use in hyperfriction later
  struct vec3d lia_v = lia_velocity(tangle, k);
  tangle->vs[k] = lia_v;

  struct vec3d evs;
  get_vs(&tangle->vnodes[k], t, &evs);
  vec3_add(&tangle->vs[k], &tangle->vs[k], &evs);

  struct octree *tree = NULL;
  if(use_BH)
    {
      //use the tree approximation to the full Biot-Savart
      tree = octree_build(tangle);
      struct vec3d v_tree;
      octree_get_vs(tree, &tangle->vnodes[k], BH_resolution, &v_tree);
    }
  else
    {
      //integrate Biot-Savart as usual
      for(m=0; m<tangle->N; ++m)
	{
	  if(tangle->connections[m].forward == -1 ||
	     m == k                               ||
	     k == tangle->connections[m].forward)
	    continue;

	  struct vec3d segment_vel = segment_field(tangle, m, tangle->vnodes[k]);
	  for(i=0; i<3; ++i)
	    tangle->vs[k].p[i] += segment_vel.p[i];
	}
    }

  //calculate the velocity due to boundary images

  struct vec3d shift_r, v_shift;
  struct vec3d v_shift_total = vec3(0, 0, 0);

  for(int j = 0; j < tangle->bimg.n; ++j) {
    shift_r = shifted(&tangle->bimg.images[j], tangle, &tangle->vnodes[k]);
    v_shift = calculate_vs(tangle, shift_r, -1);
    if(tangle->bimg.images[j].reflect > -1) //mirror wall
	    v_shift = mirror_dir_reflect(&v_shift, tangle->bimg.images[j].reflect);
    vec3_add(&v_shift_total, &v_shift_total, &v_shift);
  }
  //add everything to the result
  vec3_add(&tangle->vs[k], &tangle->vs[k], &v_shift_total);

  tangle->vels[k] = tangle->vs[k];

  if(hyperfriction) {
    //hyperfriction damps segments of vortices with very high curvature, idea is similar to hyperviscosity
    //dissipation similar to mutual friction is applied with high alpha (~0.5 -- 1)
    //but only to parts where curvature is higher than settable threshold
    //and only through the (curvature-sensitive) locally induced velocity 
    struct vec3d tmp;
    double spp = vec3_d(&tangle->normals[k]);
    if(spp > max_curvature_scale/global_dl_max) {
      vec3_cross(&tmp, &tangle->tangents[k], &lia_v);
      vec3_mul(&tmp, &tmp, -hyperalpha);
      vec3_add(&tangle->vels[k], &tangle->vels[k], &tmp);
    }
  }

  if(use_mutual_friction) {
    struct vec3d tmp, dv;

    //the velocity difference
    get_vn(&tangle->vnodes[k], t, &dv);
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

  if(tangle->status[k].status == PINNED_SLIP) {
    struct vec3d n = boundary_normals[tangle->status[k].pin_wall];
    double normal_velocity = vec3_dot(&n, &tangle->vels[k]);
    vec3_mul(&n, &n, normal_velocity);
    vec3_sub(&tangle->vels[k], &tangle->vels[k], &n);
  }
}

void update_velocities(struct tangle_state *tangle, double t)
{
  int i;
  #pragma omp parallel private(i) num_threads(global_num_threads)
  {
    #pragma omp for
    for(i=0; i<tangle->N; ++i) {
	    update_velocity(tangle, i, t);
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
     tangle->status[tangle->next_free].status == EMPTY)
    return tangle->next_free++;
  
  //otherwise we have to search for it

  for(int k=0; k<tangle->N; ++k)
    if(tangle->status[k].status == EMPTY)
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
    if(tangle->status[k].status == EMPTY)
      sum++;
  return sum;
}

static inline int out_of_box(const struct tangle_state *tangle, const struct vec3d *where)
{
  const double x = where->p[0];
  const double y = where->p[1];
  const double z = where->p[2];
  const double bounds[] = {
    tangle->box.bottom_left_back.p[0],
    tangle->box.top_right_front.p[0],
    tangle->box.bottom_left_back.p[1],
    tangle->box.top_right_front.p[1],
    tangle->box.bottom_left_back.p[2],
    tangle->box.top_right_front.p[2]
  };

  int face = -1;
  if(x < bounds[X_L])
    face = X_L;
  if(x > bounds[X_H])
    face = X_H;
  if(y < bounds[Y_L])
    face = Y_L;
  if(y > bounds[Y_H])
    face = Y_H;
  if(z < bounds[Z_L])
    face = Z_L;
  if(z > bounds[Z_H])
    face = Z_H;

  if(face > -1 && tangle->box.wall[face] != WALL_PERIODIC)
    face = -1;

  return face;
}

void enforce_boundaries(struct tangle_state *tangle)
{
  int face;
  for(int k=0; k < tangle->N; ++k)
    {
      if(tangle->status[k].status == EMPTY)
	continue;
      face = out_of_box(tangle, &tangle->vnodes[k]);
      if(face >= 0)
	  {
	    //this should be only possible with periodic faces
	    assert_msg(tangle->status[k].status != PINNED ||
		       tangle->status[k].status != PINNED_SLIP,
		       "pinned node outside of the box\n"
		       "this should have been caught with reconnections")
	    while((face=out_of_box(tangle, &tangle->vnodes[k])) >= 0)
	      tangle->vnodes[k] = periodic_shift(&tangle->vnodes[k], &tangle->box, face);
	  }
    }
}

void remove_point(struct tangle_state *tangle, int point_idx);
int add_point(struct tangle_state *tangle, int point_idx);
void remesh(struct tangle_state *tangle, double min_dist, double max_dist)
{
  int added = 0;
  for(int k=0; k<tangle->N; ++k) {
    if(tangle->status[k].status == EMPTY)
	    continue;


    int next = tangle->connections[k].forward;
    int prev = tangle->connections[k].reverse;

    struct segment sf;
    struct segment sr;
    double lf = 0;
    double lr = 0;

    if(next >= 0) {
      sf = seg_pwrap(&tangle->vnodes[k], &tangle->vnodes[next], &tangle->box);
      lf = segment_len(&sf);
    }

    if(prev >= 0) {
      sr = seg_pwrap(&tangle->vnodes[k], &tangle->vnodes[prev], &tangle->box);
      lr = segment_len(&sr);
    }

    //can we remove point k?
    if(next >= 0 && prev >= 0 && ((lf < min_dist || lr < min_dist) && (lf + lr) < max_dist )) {
      remove_point(tangle, k);
    }

    //do we need an extra point?
    if(next >= 0 && (lf > max_dist)) { //since we are adding between k and next, check only lf
      added++;
      int new_pt = add_point(tangle, k);
      update_tangent_normal(tangle, new_pt);
    }
  }
  //we could have added points outside of the domain
  if(added)
    enforce_boundaries(tangle);
}

void eliminate_loops_near_origin(struct tangle_state *tangle, double cutoff);
void eliminate_loops_near_zaxis(struct tangle_state *tangle, double cutoff);
void eliminate_small_loops(struct tangle_state *tangle, int loop_length)
{
  if(eliminate_origin_loops)
    eliminate_loops_near_origin(tangle, eliminate_loops_origin_cutoff);
  if(eliminate_zaxis_loops)
    eliminate_loops_near_zaxis(tangle, eliminate_loops_zaxis_cutoff);

  int killed = 0;

  for(int k=0; k < tangle->N; ++k)
    tangle->recalculate[k] = 0;

  for(int k=0; k < tangle->N; ++k) {
    if(tangle->status[k].status == EMPTY || tangle->recalculate[k])
	    continue; //empty or already visited point

    tangle->recalculate[k]++;

    int loop = 0;
    int here = k;
    int next = tangle->connections[here].forward;
    while(next != k) {
	    if((tangle->status[here].status == PINNED || tangle->status[here].status == PINNED_SLIP) && next < 0) {
	      //we hit a wall, turn back from k
	      here = k;
	      next = tangle->connections[here].reverse;
        while(next >= 0) {
          tangle->recalculate[here]++;
          here = next;
          next = tangle->connections[here].reverse;
          loop++;
		    }
	      //'here' now points to a node with its back to the wall
	      break;//exit the outer loop,
	    }
      tangle->recalculate[here]++;
      here = next;
      next = tangle->connections[here].forward;
      loop++;
	  } 
    if(loop < loop_length) { //the loop is short, delete it
      /*
      * for loops, the starting point doesn't matter
      * but for wall-pinned lines the code bellow only goes
      * forward, so we have to start at the end facing away from
      * the wall
      */
      killed++;
      next = here;
      while(1) {
        int tmp = next;
        next = tangle->connections[next].forward;
        printf("deleting %d\n", tmp);
        tangle->connections[tmp].forward = -1;
        tangle->connections[tmp].reverse = -1;
        tangle->status[tmp].status = EMPTY;
        if(next == here || next < 0)
          break;
      }
	  }
  }
  //printf("Killed %d loops.\n", killed);
}

void eliminate_loops_near_origin(struct tangle_state *tangle, double cutoff)
{
  for(int k=0; k < tangle->N; ++k)
    tangle->recalculate[k] = 0;

  for(int k=0; k < tangle->N; ++k)
    {
      if(tangle->status[k].status == EMPTY ||
	 tangle->recalculate[k] > 0)
	continue;

      int here = tangle->connections[k].forward;
      int next = tangle->connections[here].forward;
      int cut = vec3_d(&tangle->vnodes[k]) < cutoff;
      while(here != k)
	{
	  cut = cut && (vec3_d(&tangle->vnodes[here]) < cutoff);
	  here = next;
	  next = tangle->connections[here].forward;
	}

      here = tangle->connections[k].forward;
      next = tangle->connections[here].forward;
      if(cut) // all the points are within cutoff
	{
	  while(here != k)
	    {
	      remove_point(tangle, here);
	      here = next;
	      next = tangle->connections[here].forward;
	    }
	  remove_point(tangle, here);
	}

    }
}

void eliminate_loops_near_zaxis(struct tangle_state *tangle, double cutoff)
{
  for(int k=0; k<tangle->N; ++k)
    tangle->recalculate[k] = 0;

  for(int k=0; k<tangle->N; ++k)
    {
      if(tangle->status[k].status == EMPTY ||
	 tangle->recalculate[k] > 0)
	continue;

      int here = tangle->connections[k].forward;
      int next = tangle->connections[here].forward;
#define Z_R(v) (sqrt(v.p[0]*v.p[0] + v.p[1]*v.p[1]))
      int cut = Z_R(tangle->vnodes[k]) < cutoff;
      while(here != k)
	{
	  cut = cut && (Z_R(tangle->vnodes[here]) < cutoff);
	  here = next;
	  next = tangle->connections[here].forward;
	}

      here = tangle->connections[k].forward;
      next = tangle->connections[here].forward;
      if(cut)
	{
	  while(here != k)
	    {
	      remove_point(tangle, here);
	      here = next;
	      next = tangle->connections[here].forward;
	    }
	  remove_point(tangle, here);
	}
    }
#undef Z_R
}

void remove_point(struct tangle_state *tangle, int point_idx)
{
  int prev = tangle->connections[point_idx].reverse;
  int next = tangle->connections[point_idx].forward;

  if(prev >= 0)
    tangle->connections[prev].forward = next;
  if(next >= 0)
    tangle->connections[next].reverse = prev;

  tangle->connections[point_idx].reverse =
    tangle->connections[point_idx].forward = -1;

  tangle->status[point_idx].status = EMPTY;
}


//add a point between p and p+1 (p+1 in the sense of connections)
int add_point(struct tangle_state *tangle, int p)
{
  int next = tangle->connections[p].forward;
  int new_pt = get_tangle_next_free(tangle);
  tangle->status[new_pt].status = FREE;
  tangle->status[new_pt].pin_wall = NOT_A_FACE;

  update_tangent_normal(tangle, p);
  update_tangent_normal(tangle, next);

  struct vec3d s0 = tangle->vnodes[p];
  struct vec3d s1 = tangle->vnodes[next];
  struct segment seg = seg_pwrap(&s0, &s1, &tangle->box);
  s1 = seg.r2;

  struct vec3d s0pp = tangle->normals[p];
  struct vec3d s1pp = tangle->normals[next];

  double l = vec3_dist(&s0, &s1);

  struct vec3d a, b, new;

  struct vec3d n;
  vec3_add(&n, &s0pp, &s1pp);
  vec3_mul(&n, &n, 0.5);
  double nd = vec3_d(&n);
  if(nd < max_curvature_scale/global_dl_max && nd > 1e-5) {
    //if the curvature is too high (can happen near walls due to numerical instability) this extrapolation places the
    //point too far and breaks the curvature calculation in the subsequent steps

    //n will be identically 0 for a straight vortex
    //for both of these cases it's better to just average s0 and s1
    double R = 1/vec3_d(&n);
    double dR = R*R - l*l/4;
    //dR can become < 0 for sharp cusps
    //simplest way to deal with it is to treat the s0 and s1 as sitting
    //on opposite ends of a circle, for which dR = 0
    //this does not preserve curvature, but this is below our resolution anyway
    double delta = dR > 0 ? R - sqrt(dR) : R;

    vec3_normalize(&n);
    vec3_mul(&n, &n, -1);

    vec3_add(&a, &s0, &s1);
    vec3_mul(&a, &a, 0.5);

    vec3_mul(&b, &n, delta);

    vec3_add(&new, &a, &b);
  }
  else { //we basically have a straight vortex, just average s0 and s1
    vec3_add(&new, &s0, &s1);
    vec3_mul(&new, &new, 0.5);
  }

  tangle->vnodes[new_pt] = new;
  tangle->connections[new_pt].reverse = p;
  tangle->connections[new_pt].forward = next;
  tangle->connections[p].forward = new_pt;
  tangle->connections[next].reverse = new_pt;
  return new_pt;
}

int curvature_smoothing(struct tangle_state *tangle, double max_spp, double damping)
{
  update_tangents_normals(tangle);

  for(int i=0; i < tangle->N; ++i) {
    if(tangle->status[i].status == EMPTY ||
       tangle->status[i].status == PINNED ||
       tangle->status[i].status == PINNED_SLIP)
      continue;

    struct vec3d c = tangle->normals[i];
    double spp = vec3_d(&c);

    if(spp > max_spp)	{
      struct vec3d lia_v = lia_velocity(tangle, i);
      //vec3_normalize(&n);
      struct vec3d n;
      vec3_mul(&n, &lia_v, damping*global_dt);
      vec3_add(&tangle->vnodes[i], &tangle->vnodes[i], &n);
	  }
  }
  return 0;
}

int tangle_total_points(struct tangle_state *tangle)
{
  int total = 0;
  for(int k=0; k<tangle->N; ++k) {
    if(tangle->status[k].status != EMPTY)
	    total++;
  }
  return total;
}
