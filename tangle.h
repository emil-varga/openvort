#ifndef TANGLE_H
#define TANGLE_H

#include <stdlib.h>

#include "vec3_maths.h"

typedef enum _bc {
  OPEN,
  PIN_WALL,
  SLIP_WALL,
  PERIODIC
} boundary_conditions;

typedef enum _fc {
  X_L, X_H,
  Y_L, Y_H,
  Z_L, Z_H
} boundary_faces;

typedef enum _ns {
  EMPTY = -1, //no point here
  FREE, //ordinary vortex point
  PINNED, //pinned on the wall
  PINNED_SLIP //pinned, but can slip
} node_status;

struct neighbour_t {
  int forward;
  int reverse;
};

struct tangle_state {
  struct vec3d *vnodes;
  struct vec3d *vnodes_new;
  struct vec3d *vels;
  struct vec3d *tangents;
  struct vec3d *normals;
  
  struct neighbour_t *connections;

  node_status *status; //pinned, free, empy, etc.
  //the boundary conditions in the 6 cardinal directions
  boundary_conditions bc[6];
  
  size_t N;
  int next_free;
  int total_free;
};

void alloc_arrays(struct tangle_state *tangle, size_t n);
void expand_arrays(struct tangle_state *tangle, size_t n);
void free_arrays(struct tangle_state *tangle);
int num_free_points(struct tangle_state *tangle);

void step_nodes2(struct tangle_state *result, const struct tangle_state *tangle, double dt);
void update_tangents_normals(struct tangle_state *tangle);
void update_velocities(struct tangle_state *tangle);

/**
  @brief calculates the superfluid velocity

  Calculates and returns the superfluid velocity v_s induced by tangle at location r.
  Optionally disables one node given by skip. Use skip=-1 to not skip anything.
*/
struct vec3d calculate_vs(struct tangle_state *tangle, vec3d r, size_t skip);

void remesh(struct tangle_state *tangle, double min_dist, double max_dist);

int get_tangle_next_free(struct tangle_state *tangle);

static inline void step_nodes(struct tangle_state *tangle, double dt)
{
  step_nodes2(tangle, tangle, dt);
}

static inline void update_tangle(struct tangle_state *tangle)
{
  update_tangents_normals(tangle);
  update_velocities(tangle);
}


#endif //TANGLE_H
