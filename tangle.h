#ifndef TANGLE_H
#define TANGLE_H

#include <stdlib.h>

#include "vec3_maths.h"

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
  
  size_t N;
  int next_free;
  int total_free;
};

void alloc_arrays(struct tangle_state *tangle, size_t n);
void free_arrays(struct tangle_state *tangle);
int num_free_points(struct tangle_state *tangle);

void step_nodes(struct tangle_state *tangle, double dt)
{
  step_nodes2(tangle, tangle, dt);
}
void step_nodes2(struct tangle_state *result, const struct tangle_state *tangle, double dt);
void update_tangents_normals(struct tangle_state *tangle);
void update_velocities(struct tangle_state *tangle);

void update_tangle(struct tangle_state *tangle)
{
  update_tangents_normals(tangle);
  update_velocities(tangle);
}

int get_tangle_next_free(struct tangle_state *tangle);

#endif //TANGLE_H
