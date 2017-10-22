#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <printf.h>
#include <stddef.h>
#include <math.h>

#include <fenv.h>

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
  feenableexcept(FE_OVERFLOW | FE_UNDERFLOW | FE_INVALID | FE_DIVBYZERO);
  struct tangle_state *tangle = (struct tangle_state*)malloc(sizeof(struct tangle_state));

  alloc_arrays(tangle, 512);
  struct vec3d center1 = vec3(0, 0, 0);
  struct vec3d dir1    = vec3(0, 0, 1);
  size_t k;
  int recs = 0;
  char filename[128];
  
  add_circle(tangle, &center1, &dir1, 0.1, 64);
  save_tangle("v1.dat", tangle);
  add_circle(tangle, &center1, &dir1, 0.09, 64);
  save_tangle("v2.dat", tangle);

  const int Nshot = 100;
  int shot = Nshot;
  int frame = 0;

  for(k=0; recs==0 && k < 100000; ++k)
    {
      printf("Step %04zu\n", k);
      update_tangle(tangle);
      if(!shot)
	{
	  sprintf(filename, "data/frame%04zu.dat", frame);
	  save_tangle(filename, tangle);
	  frame++;
	  shot = Nshot;
	} 
      //      if(reconnect(tangle, 1e-3, 0) > 0)
      //printf("reconnected!\n");
      rk4_step(tangle, 1e-3);
      shot--;
    }
  free_arrays(tangle);
  free(tangle);
}
