#include <stdio.h>
#include <stdlib.h>
#include <printf.h>
#include <stddef.h>
#include <math.h>

#define _DEBUG_

#include "vec3_maths.h"
#include "tangle.h"
#include "vortex_utils.h"

void vec3_print(const struct vec3d *v)
{
  printf("(%g, %g, %g)", v->p[0], v->p[1], v->p[2]);
}

int main(int argc, char **argv)
{
  struct tangle_state *tangle = (struct tangle_state*)malloc(sizeof(struct tangle_state));

  alloc_arrays(tangle, 256);
  struct vec3d center = vec3(0, 0, 0);
  struct vec3d dir    = vec3(0, 1, 0);
  add_circle(tangle, &center, &dir, 1, 128);
  update_tangle(tangle);
  save_tangle("steapa.dat", tangle);
  remesh(tangle, 0, 3e-2);
  printf("saving step0...");
  save_tangle("step0.dat", tangle);
  printf("done\n");
  update_tangents_normals(tangle);
  update_velocities(tangle);
  step_nodes(tangle, 1e-3);
  save_tangle("step1.dat", tangle);
  free_arrays(tangle);
  free(tangle);
}
