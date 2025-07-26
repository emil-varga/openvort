/*******************************************************************************
 * Copyright (C) 2018 Emil Varga <varga.emil@gmail.com>
 * 
 * This file is part of OpenVort
 * 
 * OpenVort is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
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
#include "vortex_injection.h"


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

  printf("Reading config from %s.\n", conf_file);
  char filename[128];
  if(!load_conf(conf_file, tangle)) {
      printf("Can't initialize! Exiting.\n");
      free_tangle(tangle);
      free(tangle);
      return EXIT_FAILURE;
  }

  print_config(tangle);

  enforce_boundaries(tangle);

  remesh(tangle, global_dl_min, global_dl_max);

  eliminate_small_loops(tangle, global_small_loop_cutoff);

  update_tangle(tangle, 0);

  int shot = global_frame_shot - 1;
  int frame = 0;
  int recs = 0;
  int nrec = 0;

  struct timespec t0, ti;
  clockid_t clock = CLOCK_MONOTONIC;
  clock_gettime(clock, &t0);
  int Np = tangle_total_points(tangle);
  double time = 0;
  fflush(stdout);
  sprintf(filename, "%s/init.dat", output_dir);
  save_tangle(filename, tangle);
  for(int k=0; Np > 0; ++k) {
    printf("Step %d, time = %g, recs: %d, Np: %d\n", k, time, recs, Np);
    nrec = reconnect(tangle, time, rec_dist, global_reconnection_angle_cutoff);

    recs += nrec;
    eliminate_small_loops(tangle, global_small_loop_cutoff);
    if(shot <= 0)	{
      update_tangle(tangle, time);
      //output_dir declared extern in configuration.h
      sprintf(filename, "%s/frame%08d.dat", output_dir, frame);
      save_tangle(filename, tangle);
      frame++;
      shot = global_frame_shot;
    }

    inject_vortices(tangle, time);
    update_tangle(tangle, time);
    rk4_step(tangle, time, global_dt);
    update_tangents_normals(tangle);
    remesh(tangle, global_dl_min, global_dl_max);
    eliminate_small_loops(tangle, global_small_loop_cutoff);
    enforce_boundaries(tangle);
    //curvature_smoothing(tangle, 0.5/global_dl_max, 0.1);

    // #ifdef _DEBUG_
    // printf("Checking tangle integrity.");
    // check_integrity(tangle);
    // #endif

    Np = tangle_total_points(tangle);
    shot--;
    time += global_dt;
    fflush(stdout);
    if(global_max_steps > 0 && k > global_max_steps)
      break;
  }
  clock_gettime(clock, &ti);
  printf("Elapsed seconds: %f\n", time_diff(&t0, &ti));
  free_tangle(tangle);
  free(tangle);

  return EXIT_SUCCESS;
}
