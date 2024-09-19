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

#include "vortex_dynamics.h"
#include "vec3_maths.h"
#include "util.h"
#include "vortex_constants.h"
#include "vortex_utils.h"
#include "tangle.h"
#include <math.h> //for cos()
#include <stdint.h>
#include <stdio.h>

/*
 * Euler-steps the tangle in 'tangle' with time step dt and saves the result to 'result'.
 *
 * Optionally can use externally supplied velocities use_vel. Pass NULL to use standard vortex
 * velocity.
 */
void euler_step2(struct tangle_state *result,
		 const struct tangle_state *tangle, double dt,
		 const struct vec3d *use_vel)
{
  int k;

  for(k=0; k<tangle->N; ++k) {
    if(tangle->status[k].status == EMPTY)
	    continue; //empty node
    struct vec3d move;
    if(!use_vel) {
      vec3_mul(&move, &tangle->vels[k], dt);
      vec3_add(&result->vnodes[k], &tangle->vnodes[k], &move);
	  } else {
	    vec3_mul(&move, &use_vel[k], dt);
	    vec3_add(&result->vnodes[k], &tangle->vnodes[k], &move);
	  }
  }
}

void euler_step(struct tangle_state *tangle, double dt)
{
  euler_step2(tangle, tangle, dt, NULL);
}

void rk4_step2(struct tangle_state *result, const struct tangle_state *tangle, double t, double dt)
{
  int N;
  struct tangle_state rk_state[3]; //k2 through k4, we already have k1 in tangle

  N = tangle->N;
  for(int k=0; k < 3; ++k) {
    create_tangle(&rk_state[k], N);
    //only copy the stuff we actualy need
    rk_state[k].bimg = tangle->bimg;
    rk_state[k].box = tangle->box;
    for(int kk=0; kk < N; ++kk) {
      rk_state[k].status[kk] = tangle->status[kk];
      rk_state[k].connections[kk] = tangle->connections[kk];
	  }
  }

  //calculate k2
  euler_step2(&rk_state[0], tangle, dt/2, NULL);
  // enforce_boundaries(&rk_state[0]);
  update_tangents_normals(&rk_state[0]);
  update_velocities(&rk_state[0], t + dt/2, NULL);

  //calculate k3
  euler_step2(&rk_state[1], tangle, dt/2, rk_state[0].vels);
  // enforce_boundaries(&rk_state[1]);
  update_tangents_normals(&rk_state[1]);
  update_velocities(&rk_state[1], t+dt/2, NULL);

  //calculate k4
  euler_step2(&rk_state[2], tangle, dt/2, rk_state[1].vels);
  // enforce_boundaries(&rk_state[2]);
  update_tangents_normals(&rk_state[2]);
  update_velocities(&rk_state[2], t + dt, NULL);

  for(int k=0; k < N; ++k) {
    if(tangle->status[k].status == EMPTY)
	    continue;

    //move contains the full step, a are partial steps
    struct vec3d move, a;
    vec3_mul(&move, &tangle->vels[k], dt/6); //k1

    vec3_mul(&a, &rk_state[0].vels[k], dt/3); // k2
    vec3_add(&move, &move, &a);

    vec3_mul(&a, &rk_state[1].vels[k], dt/3); // k3
    vec3_add(&move, &move, &a);

    vec3_mul(&a, &rk_state[2].vels[k], dt/6); //k4
    vec3_add(&move, &move, &a);

    vec3_add(&result->vnodes[k], &result->vnodes[k], &move);
  }

  for(int k=0; k<3; ++k) {
    free_tangle(&rk_state[k]);
  }
}

void rk4_step(struct tangle_state *tangle, double t, double dt)
{
  rk4_step2(tangle, tangle, t, dt);
}

//helper function to swap around the connection indices
int do_reconnection(struct tangle_state *tangle, size_t k, size_t l);

//helper for wall-reconnections
int check_wall(struct tangle_state *tangle, int k, int wall, double rdist);
int connect_to_wall(struct tangle_state *tangle, int k, int wall, node_status_t pin_mode, double time);

/*
  Run through all the pairs of nodes and check their distance if they are not
  immediate neighbors. Reconnect them if they are close enough. Also reconnect points close
  to a solid wall.

  This could, and should, in the future include some smarter estimation of possibility
  of a reconnections -- i.e., where in the BSP tree the two nodes are and only check them
  if can, in principle, be close enough.
 */
int reconnect(struct tangle_state *tangle, double t, double rec_dist, double rec_angle)
{
  int k, l;
  struct vec3d *v1, *v2; //points under test
  struct vec3d d1, d2; //direction vectors from v1, v2
  double calpha; //cosine of alpha between d1, d2
  int Nrecs = 0;

  for(k=0; k<tangle->N; ++k)
    tangle->recalculate[k] = 0;

  /*
   * "reconnect" with walls
   */
  for(int wall = 0; wall < 6; ++wall) {
    if(tangle->box.wall[wall] == WALL_OPEN || tangle->box.wall[wall] == WALL_PERIODIC)
	    continue;

    for(k=0; k<tangle->N; ++k) {
	    if(tangle->status[k].status == PINNED       ||
	        tangle->status[k].status == PINNED_SLIP ||
	        tangle->status[k].status == EMPTY       ||
	        tangle->recalculate[k] > 0)
	      continue;

	    if(check_wall(tangle, k, wall, rec_dist/2)) {
        Nrecs += connect_to_wall(tangle, k, wall, pin_mode, t);
	    }
	  }
  }

  /*
   * Wall reconnection can create small loops that could be broken by the domain killing below
   * which can cause the simulation to crash. These loops need to be removed.
   */
  if(Nrecs > 0) {
      eliminate_small_loops(tangle, small_loop_cutoff);
      //update_tangents_normals(tangle);
  }

  /*
   * Eliminate all points that are active and outside the walls. This can potentially happen
   * if the time step is too long and more than one discretisation point gets outside the domain.
   */
  int domain_killed = 0;
  for(int wall = 0; wall < 6; ++wall) {
    if(tangle->box.wall[wall] == WALL_OPEN || tangle->box.wall[wall] == WALL_PERIODIC)
	    continue;

    for(k=0; k<tangle->N; ++k) {
	    if(tangle->status[k].status == EMPTY)
	      continue;

      if(wall_dist(tangle, k, wall) < 0) {
          remove_point(tangle, k, 0);
          domain_killed++;
      }
	  }
  }

   /*
   * Some small loops might have been created by the removals above.
   * These can confuse the tangent calculation so they should be removed.
   */
  if(domain_killed > 0) {
    printf("Removed %d points outside the domain.\n", domain_killed);
    eliminate_small_loops(tangle, small_loop_cutoff);
    //update_tangents_normals(tangle);
  }

  /*
   * Now do standard vortex-vortex reconnections
   */
  //Nrecs = 0;
  for(k=0; k < tangle->N; ++k)
    tangle->recalculate[k] = 0;

  for(k=0; k<tangle->N; ++k) {
    //do not reconnect points attached to the walls

    if(tangle->status[k].status == EMPTY       ||
       tangle->status[k].status == PINNED      ||
	     tangle->status[k].status == PINNED_SLIP ||
	     tangle->recalculate[k])
	    continue; //skip empty nodes and nodes that reconnected in this pass
      
    for(l=k+1; l<tangle->N; ++l) {
	    if(tangle->status[l].status == EMPTY       ||
	       tangle->status[l].status == PINNED      ||
	       tangle->status[l].status == PINNED_SLIP ||
	       tangle->connections[k].forward == l     ||
	       tangle->connections[k].reverse == l     ||
	       tangle->recalculate[l])
	      continue; //skip empty nodes and neighbors of k
	                //and nodes that went through a reconnection 
	    struct segment seg = seg_pwrap(tangle->vnodes + k, tangle->vnodes + l, &tangle->box);
	    v1 = &seg.r1;
	    v2 = &seg.r2;

	    if(vec3_dist(v1, v2) > rec_dist)
	      continue;
	  
      //the nodes are close and they are not neighbors
      //now check the angle

      //update_tangents_normals(tangle);

      d1 = tangle->tangents[k];
      d2 = tangle->tangents[l];

      //In some cases vortex-vortex reconnections can create tiny loops (n ~ 3)
      //pinned on the wall which have colinear points for which calculation of the tangents fails.
      //These result in tangents of 0 length. If the tangent (which should be 1) is small
      //we just ignore this point because this loop will be removed on the next sweep of eliminate_small_loops
      //The 0.5 threshold is arbitrary. It will always be roughly 1 or roughly 0.
      if(vec3_d(&d1) < 0.5 || vec3_d(&d2) < 0.5)
        continue;
      //normalized dot -- just the cosine of the angle between vectors
      calpha = vec3_ndot(&d1, &d2);

      //We want to reconnect for angles LARGER than rec_angle
      //(i.e., parallel lines do not reconnect) and cosine is
      //a decreasing function.
      if(calpha > cos(rec_angle))
        continue; //angle too small

      //Do not reconnect if the nodes are getting further apart from each other
      struct vec3d dx, dv;
      //we only update velocities if we really need them
      //do not bother with the Barnes-Hut here, passing NULL for the tree turns it off even if enabled globally
      update_velocity(tangle, k, t, NULL);
      update_velocity(tangle, l, t, NULL);
      vec3_sub(&dx, v1, v2);
      vec3_sub(&dv, &tangle->vels[k], &tangle->vels[l]);

      if(vec3_dot(&dx, &dv) > 0)
        continue; //the points are getting further from each other

      //save the old connections because we might need them later
      int knext = tangle->connections[k].forward;
      int kprev = tangle->connections[k].reverse;
      int lnext = tangle->connections[l].forward;
      int lprev = tangle->connections[l].reverse;

      //angle is not too close and we can finally reconnect points k and l
      if(!do_reconnection(tangle, k, l))
        continue; //do_reconnect additionally checks if the total length doesn't increase

      //flag the neighbourhood as tainted so that it doesn't flash
      //back and forth in a single pass
      tangle->recalculate[k]++;
      tangle->recalculate[tangle->connections[k].forward]++;
      tangle->recalculate[tangle->connections[k].reverse]++;
      tangle->recalculate[knext]++;
      tangle->recalculate[kprev]++;

      tangle->recalculate[l]++;
      tangle->recalculate[tangle->connections[l].forward]++;
      tangle->recalculate[tangle->connections[l].reverse]++;
      tangle->recalculate[lnext]++;
      tangle->recalculate[lprev]++;

      Nrecs++;
      //printf("vortex-vortex reconnection %d %d\n", k, l);
      
      break; 
	  }
  }
  return Nrecs;
}

int do_reconnection(struct tangle_state *tangle, size_t k, size_t l)
{
  size_t kf, kr, lf, lr;

  double cf1, cf2;

  kf = tangle->connections[k].forward;
  kr = tangle->connections[k].reverse;

  lf = tangle->connections[l].forward;
  lr = tangle->connections[l].reverse;

  #define PSEG(x,y) seg_pwrap(tangle->vnodes+x, tangle->vnodes+y, &tangle->box)
  struct segment seg_kl   = PSEG(k, l);
  struct segment seg_kflr = PSEG(kf, lr);
  struct segment seg_kkf  = PSEG(k, kf);
  struct segment seg_llr  = PSEG(l, lr);

  struct segment seg_krlf = PSEG(kr, lf);
  struct segment seg_kkr  = PSEG(k, kr);
  struct segment seg_llf  = PSEG(l, lf);
  #undef PSEG

  cf1 = segment_len(&seg_kl) +
      segment_len(&seg_kflr) -
      segment_len(&seg_kkf)  -
      segment_len(&seg_llr);

  cf2 = segment_len(&seg_kl) +
      segment_len(&seg_krlf) -
      segment_len(&seg_kkr)  -
      segment_len(&seg_llf);

  if(cf1 > 0 && cf2 > 0)
    return 0; //the reconnection increases the size, don't do it;

  if(cf1 < cf2)
    {
      tangle->connections[k].forward = l;
      tangle->connections[l].reverse = k;

      tangle->connections[kf].reverse = lr;
      tangle->connections[lr].forward = kf;
    }
  else
    {
      tangle->connections[l].forward = k;
      tangle->connections[k].reverse = l;

      tangle->connections[lf].reverse = kr;
      tangle->connections[kr].forward = lf;
    }

  return 1;
}

int check_wall(struct tangle_state *tangle, int k, int wall, double rdist)
{
  /*
   * Checks whether a node k is closer than rdist to the wall.
   */
  int idx[6];
  idx[X_L] = idx[X_H] = 0;
  idx[Y_L] = idx[Y_H] = 1;
  idx[Z_L] = idx[Z_H] = 2;
  switch(wall)
  {
    case X_L:
    case Y_L:
    case Z_L:
      return tangle->vnodes[k].p[idx[wall]] < (rdist + tangle->box.bottom_left_back.p[idx[wall]]);

    case X_H:
    case Y_H:
    case Z_H:
      return tangle->vnodes[k].p[idx[wall]] > (tangle->box.top_right_front.p[idx[wall]] - rdist);

    default:
      break;
  }
  return 0;
}

int connect_to_wall(struct tangle_state *tangle, int k, int wall, node_status_t pin_mode, double time)
{
  update_velocity(tangle, k, time, NULL); //do not bother with B-H for a small number of points
  double d0 = wall_dist(tangle, k, wall);
  //check that the node is actually getting closer to the wall
  //boundary_normals are inward-facing unit normals
  if(d0 > 0 && vec3_dot(&tangle->vels[k], &boundary_normals[wall]) > 0)
    return 0; //EXIT because we are getting further from the wall

  int next = tangle->connections[k].forward;
  int prev = tangle->connections[k].reverse;
  double d1 = wall_dist(tangle, next, wall);
  double dm1 = wall_dist(tangle, prev, wall);

  struct segment s1 = seg_pwrap(&tangle->vnodes[k], &tangle->vnodes[next], &tangle->box);
  struct segment sm1 = seg_pwrap(&tangle->vnodes[k], &tangle->vnodes[prev], &tangle->box);
  double vd1 = segment_len(&s1);
  double vdm1 = segment_len(&sm1);

  //check that the total length does not increase
  if(vd1 < d0 + d1 && vdm1 < d0 + dm1)
    return 0;

  //printf("Pinning %d %d %d\n", prev, k, next);

  struct vec3d tmp;

  //point 1
  int new_pt = get_tangle_next_free(tangle);
  tangle->status[new_pt].status = pin_mode;
  tangle->status[new_pt].pin_wall = wall;

  //point 2
  int new_pt2 = get_tangle_next_free(tangle);
  tangle->status[new_pt2].status = pin_mode;
  tangle->status[new_pt2].pin_wall = wall;
  //printf("New points: %d %d\n", new_pt, new_pt2);

  //the first new point is the projection of the kth point on the wall
  vec3_mul(&tmp, &boundary_normals[wall], -d0);
  vec3_add(&tangle->vnodes[new_pt], &tangle->vnodes[k], &tmp);

  //flag the changed points to avoid overcrowding reconnections
  tangle->recalculate[k]++;
  tangle->recalculate[new_pt]++;
  tangle->recalculate[new_pt2]++;
  tangle->recalculate[next]++;
  tangle->recalculate[prev]++;

  if(d1 < dm1) {
    vec3_mul(&tmp, &boundary_normals[wall], -d1);
    vec3_add(&tangle->vnodes[new_pt2], &tangle->vnodes[next], &tmp);

    tangle->connections[k].forward = new_pt;
    tangle->connections[new_pt].forward = -1;
    tangle->connections[new_pt].reverse = k;
    tangle->connections[new_pt2].forward = next;
    tangle->connections[new_pt2].reverse = -1;
    tangle->connections[next].reverse = new_pt2;
  } else {
    vec3_mul(&tmp, &boundary_normals[wall], -dm1);
    vec3_add(&tangle->vnodes[new_pt2], &tangle->vnodes[prev], &tmp);

    tangle->connections[k].reverse = new_pt;
    tangle->connections[new_pt].forward = k;
    tangle->connections[new_pt].reverse = -1;
    tangle->connections[new_pt2].forward = -1;
    tangle->connections[new_pt2].reverse = prev;
    tangle->connections[prev].forward = new_pt2;
  }
  //printf("vortex-wall reconnection %d (%d, %d)\n", k, prev, next);
  return 2;
}
