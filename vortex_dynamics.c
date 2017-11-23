#include "vortex_dynamics.h"
#include "vec3_maths.h"
#include <math.h> //for cos()
#include <stdint.h>
#include <stdio.h>

void euler_step2(struct tangle_state *result,
		 const struct tangle_state *tangle, double dt)
{
  size_t k;

  for(k=0; k<tangle->N; ++k)
    {
      if(tangle->connections[k].forward==-1)
	continue; //empty node
      struct vec3d move;
      vec3_mul(&move, &tangle->vels[k], dt);
      vec3_add(&result->vnodes[k], &tangle->vnodes[k],
	       &move);
    }
}

void euler_step(struct tangle_state *tangle, double dt)
{
  euler_step2(tangle, tangle, dt);
}

void rk4_step2(struct tangle_state *result,
	       const struct tangle_state *tangle, double dt)
{
  size_t N;
  struct tangle_state rk_state[3]; //k2 through k4, we already have k1 in tangle

  N = tangle->N;
  for(int k=0; k < 3; ++k)
    alloc_arrays(&rk_state[k], N);

  //calculate k1
  euler_step2(&rk_state[0], tangle, dt/2);
  update_tangents_normals(&rk_state[0]);
  update_velocities(&rk_state[0]);

  //calculate k2
  euler_step2(&rk_state[1], &rk_state[0], dt/2);
  update_tangents_normals(&rk_state[1]);
  update_velocities(&rk_state[1]);

  //calculate k3
  euler_step2(&rk_state[2], &rk_state[1], dt/2);
  update_tangents_normals(&rk_state[2]);
  update_velocities(&rk_state[2]);

  for(int k=0; k < N; ++k)
    {
      if(tangle->connections[k].forward == -1)
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

      vec3_add(&result->vnodes[k], &result->vnodes[k],
	       &move);
    }

  for(int k=0; k<3; ++k)
    free_arrays(&rk_state[k]);
}

void rk4_step(struct tangle_state *tangle, double dt)
{
  rk4_step2(tangle, tangle, dt);
}

//helper function to swap around the connection indices
void do_reconnection(struct tangle_state *tangle, size_t k, size_t l);
size_t reconnect(struct tangle_state *tangle, double rec_dist, double rec_angle)
{
  /*
    Run through all the pairs of nodes and check their distence if they are not
    immediate neighbours. Reconnect them if they are close enough.

    This could, and should, in the future include some smarter estimation of possibility
    of a reconnections -- i.e., where in the BSP tree the two nodes are and only check them
    if can, in principle, be close enough.
   */
  size_t k, l;
  struct vec3d *v1, *v2; //points under test
  struct vec3d d1, d2; //direction vectors from v1, v2
  double calpha; //cosine of alpha between d1, d2
  size_t Nrecs = 0;

  for(k=0; k<tangle->N; ++k)
    {
      if(tangle->connections[k].forward == -1)
	continue; //skip empty nodes
      
      for(l=k+1; l<tangle->N; ++l)
	{
	  if(tangle->connections[l].forward == -1 ||
	     tangle->connections[k].forward == l  ||
	     tangle->connections[k].reverse == l)
	    continue; //skip empty nodes and neighbours of k

	  v1 = &tangle->vnodes[k];
	  v2 = &tangle->vnodes[l];

	  if(vec3_dist(v1, v2) > rec_dist)
	    continue;
	  //the nodes are close and they are not neighbouds
	  //now check the angle

	  //calculate the direction vectors d1, d2
	  vec3_sub(&d1, &tangle->vnodes[tangle->connections[k].forward],
		   v1);
	  vec3_sub(&d2, &tangle->vnodes[tangle->connections[l].forward],
		   v2);
	  //normalized dot -- just the cosine of the angle between vectors
	  calpha = vec3_ndot(&d1, &d2);

	  //We want to reconnect for angles LARGER than rec_angle
	  //(i.e., parallel lines do not reconnect) and cosine is
	  //a decreasing function.
	  if(calpha > cos(rec_angle))
	    continue; //angle too close

	  //angle is not too close and we can finally reconnect points k and l
	  do_reconnection(tangle, k, l);
	  //printf("reconnecting %d %d\n", k, l);
	  Nrecs++;
	  
	  break; //TODO: this is to allow for only one reconnection on point k
	         //this kind of solution is rather arbitrary
	}
    }
  return Nrecs;
}

void do_reconnection(struct tangle_state *tangle, size_t k, size_t l)
{
  size_t kf, lr;

#ifdef _DEBUG_
  printf("reconnecting %zu %zu\n", k, l);
#endif

  kf = tangle->connections[k].forward;
  lr = tangle->connections[l].reverse;

  tangle->connections[k].forward = l;
  tangle->connections[l].reverse = k;

  tangle->connections[lr].forward = kf;
  tangle->connections[kf].reverse = lr;
}
