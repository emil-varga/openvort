#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <printf.h>
#include <stddef.h>
#include <math.h>
#include <time.h>

#include <omp.h>
#include <fenv.h>

#include "vec3_maths.h"
#include "tangle.h"
#include "vortex_utils.h"
#include "vortex_dynamics.h"
#include "vortex_constants.h"
#include "util.h"
#include "configuration.h"


void vec3_print(const struct vec3d *v)
{
  printf("(%g, %g, %g)", v->p[0], v->p[1], v->p[2]);
}

double time_diff(struct timespec *t0, struct timespec *t1)
{
  double ns0 = t0->tv_sec*1e9 + t0->tv_nsec;
  double ns1 = t1->tv_sec*1e9 + t1->tv_nsec;

  return (ns1 - ns0)/1e9;
}

int main(int argc, char **argv)
{
  //this populates char conf_file[] and char output_dir[]
  if(!parse_options(argc, argv))
    return EXIT_FAILURE;
  setup_outdir(output_dir);

  srand48(time(NULL));
  feenableexcept(FE_OVERFLOW | FE_UNDERFLOW | FE_INVALID | FE_DIVBYZERO);
  struct tangle_state *tangle = (struct tangle_state*)malloc(sizeof(struct tangle_state));

  create_tangle(tangle, 512);

  char filename[128];
  if(!load_conf(conf_file, tangle))
    {
      free_tangle(tangle);
      free(tangle);
      return EXIT_FAILURE;
    }

  save_tangle("v0.dat", tangle);
  enforce_boundaries(tangle);
  save_tangle("v1.dat", tangle);

  const int Nshot = 10;
  int shot = Nshot - 1;
  int frame = 0;
  int recs = 0;
  int nrec = 0;

  struct timespec t0, ti;
  clockid_t clock = CLOCK_MONOTONIC;

  clock_gettime(clock, &t0);
  int Np = tangle_total_points(tangle);
  for(int k=0; Np > 0; ++k)
    {
      printf("Step %d, recs: %d, Np: %d\n", k, recs, Np);
      update_tangle(tangle);
      check_integrity(tangle);
      nrec = reconnect(tangle, 2.5e-3, DEG2RAD(5));
      check_integrity(tangle);
      eliminate_small_loops(tangle, 5);
      check_integrity(tangle);
      recs += nrec;

      if(!shot)
	{
	  //output_dir declared extern in configuration.h
	  sprintf(filename, "%s/frame%04d.dat", output_dir, frame);
	  save_tangle(filename, tangle);
	  frame++;
	  shot = Nshot;
	} 
      rk4_step(tangle, 5e-5);
      remesh(tangle, 1e-3, 5e-3);
      enforce_boundaries(tangle);

      check_integrity(tangle);
      Np = tangle_total_points(tangle);
      shot--;
    }
  clock_gettime(clock, &ti);
  printf("Elapsed seconds: %f\n", time_diff(&t0, &ti));
  free_tangle(tangle);
  free(tangle);

  return EXIT_SUCCESS;
}
