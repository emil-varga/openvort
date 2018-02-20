#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <printf.h>
#include <stddef.h>
#include <math.h>
#include <time.h> //for clock_gettime

#include <libconfig.h>

#include <omp.h>
#include <fenv.h>

#include "vec3_maths.h"
#include "tangle.h"
#include "vortex_utils.h"
#include "vortex_dynamics.h"
#include "vortex_constants.h"
#include "util.h"

#define DEG2RAD(X) ((X)*M_PI/180.0)

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
  if(argc != 2)
    {
      print_usage(argv[0]);
      return EXIT_FAILURE;
    }

  config_t cfg;
  config_setting_t *setting;

  config_init(&cfg);
  if(!config_read_file(&cfg, argv[1]))
    {
      fprintf(stderr, "Can't read config file %s", argv[1]);
      return EXIT_FAILURE;
    }

  feenableexcept(FE_OVERFLOW | FE_UNDERFLOW | FE_INVALID | FE_DIVBYZERO);
  struct tangle_state *tangle = (struct tangle_state*)malloc(sizeof(struct tangle_state));

  create_tangle(tangle, 512);
  for(int k=0; k<6; ++k)
      tangle->box.wall[k] = WALL_PERIODIC;

  struct vec3d center1 = vec3(0, 0, 0);
  struct vec3d center2 = vec3(0.05, 0, -0.01);
  struct vec3d dir1    = vec3(0, 0, 1);
  struct vec3d dir2    = vec3(0, 0, -1);
  char filename[128];

  add_circle(tangle, &center1, &dir1, 0.05, 128);
  save_tangle("v1.dat", tangle);
  enforce_boundaries(tangle);

  add_circle(tangle, &center2, &dir2, 0.09, 128);
  save_tangle("v2.dat", tangle);

  const int Nshot = 100;
  int shot = Nshot - 1;
  int frame = 0;
  int recs = 0;
  int nrec = 0;

  struct timespec t0, ti;
  clockid_t clock = CLOCK_MONOTONIC;

  clock_gettime(clock, &t0);
  for(int k=0; k < 10000; ++k)
    {
      printf("Step %d, recs: %d\n", k, recs);
      nrec = reconnect(tangle, 2.5e-3, DEG2RAD(5));
      check_integrity(tangle);
      eliminate_small_loops(tangle, 5);
      check_integrity(tangle);
      recs += nrec;

      update_tangle(tangle);
      check_integrity(tangle);
      if(!shot)
	{
	  sprintf(filename, "data_rec/frame%04d.dat", frame);
	  save_tangle(filename, tangle);
	  frame++;
	  shot = Nshot;
	} 
      rk4_step(tangle, 1e-3);
      remesh(tangle, 2.5e-3, 10e-3);
      enforce_boundaries(tangle);

      check_integrity(tangle);
      shot--;
    }
  clock_gettime(clock, &ti);
  printf("Elapsed seconds: %f\n", time_diff(&t0, &ti));
  free_tangle(tangle);
  free(tangle);
  config_destroy(&cfg);

  return EXIT_SUCCESS;
}
