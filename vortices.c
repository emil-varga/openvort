#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <printf.h>
#include <stddef.h>
#include <math.h>

#include <fenv.h>

#include "vec3_maths.h"
#include "tangle.h"
#include "vortex_utils.h"
#include "vortex_dynamics.h"
#include "vortex_constants.h"

void vec3_print(const struct vec3d *v)
{
  printf("(%g, %g, %g)", v->p[0], v->p[1], v->p[2]);
}

int main(int argc, char **argv)
{
  feenableexcept(FE_OVERFLOW | FE_UNDERFLOW | FE_INVALID | FE_DIVBYZERO);
  struct tangle_state *tangle = (struct tangle_state*)malloc(sizeof(struct tangle_state));

  alloc_arrays(tangle, 512);
  struct vec3d center1 = vec3(0, 0, 0);
  struct vec3d center2 = vec3(0.05, 0, -0.05);
  struct vec3d dir1    = vec3(0, 0, 1);
  struct vec3d dir2    = vec3(0, 0, -1);
  char filename[128];
  
  add_circle(tangle, &center1, &dir1, 0.1, 64);
  save_tangle("v1.dat", tangle);
  add_circle(tangle, &center2, &dir2, 0.09, 64);
  save_tangle("v2.dat", tangle);

  int Nshot = 100;
  int shot = Nshot;
  int frame = 0;
  int recs = 0;

  for(int k=0; k < 100000; ++k)
    {
      printf("Step %d, recs: %d\n", k, recs);
      recs += reconnect(tangle, 1.5e-2, 0.1);
      if (recs > 0) {
	Nshot = 0;
	shot = 0;
      }
      update_tangle(tangle);
      if(!shot)
	{
	  sprintf(filename, "data_rec/frame%04zu.dat", frame);
	  save_tangle(filename, tangle);
	  frame++;
	  shot = Nshot;
	} 
      rk4_step(tangle, 1e-3);
      shot--;
    }
  free_arrays(tangle);
  free(tangle);
}
