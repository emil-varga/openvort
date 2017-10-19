#ifndef VORTEX_DYNAMICS
#define VORTEX_DYNAMICS

#include "tangle.h"

/*
  All functions dealing with movement of the vortices should go here.
  
  Particularly, stepping and reconnecting.
 */


// Stepping
struct rk4_state {
  //the tangle states here should be reduced
  //only calculate what we really need for the velocity
  struct tangle_state *s1, *s2, *s3, *s4;

  double dt;
};

struct ab4_state {
  struct tangle_state *s1, *s2, *s3, *s4;

  double dt;
};

void euler_step(struct tangle_state *tangle, double dt);

void rk4_step(struct tangle_state *tangle, double dt);
void ab4_init_state(struct ab4_state **state,
		    const struct tangle_state *tangle, double dt);
void ab4_step(struct tangle_state *tangle, struct ab4_state *state);


//reconnections

size_t reconnect(struct tangle_state *tangle, double rec_dist, double rec_angle);

#endif //VORTEX_DYNAMICS
