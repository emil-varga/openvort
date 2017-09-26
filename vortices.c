#include <stdio.h>
#include <stdlib.h>
#include <printf.h>
#include <stddef.h>
#include <math.h>

#define _DEBUG_

#include "vec3_maths.h"
#include "tangle.h"
#include "vortex_utils.h"
#include "vortex_dynamics.h"

void vec3_print(const struct vec3d *v)
{
  printf("(%g, %g, %g)", v->p[0], v->p[1], v->p[2]);
}

int main(int argc, char **argv)
{
  struct tangle_state *tangle = (struct tangle_state*)malloc(sizeof(struct tangle_state));

  alloc_arrays(tangle, 512);
  struct vec3d center1 = vec3(0, 0, 0);
  struct vec3d dir1    = vec3(0, 0, 1);
  struct vec3d center2 = vec3(0.5, 0, 0.03);
  struct vec3d dir2    = vec3(0, 0, -1);
  size_t k;
  int recs = 0;
  char filename[128];
  
  add_circle(tangle, &center1, &dir1, 0.5, 256);
  save_tangle("v1.dat", tangle);
  add_circle(tangle, &center2, &dir2, 0.5, 256);
  save_tangle("v2.dat", tangle);

  for(k=0; recs==0 && k < 510; ++k)
    {
      sprintf(filename, "data/step%04zu.dat", k);
      printf("Step %04zu\n", k);
      save_tangle(filename, tangle);
      if(reconnect(tangle, 1e-3, 0) > 0)
	printf("reconnected!\n");
      update_tangle(tangle);
      step_nodes(tangle, 5e-6);
    }
  free_arrays(tangle);
  free(tangle);
}
