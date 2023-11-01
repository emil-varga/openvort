/*******************************************************************************
 * Copyright (C) 2019 Emil Varga <varga.emil@gmail.com>
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


#include "tangle.h"
#include "vortex_utils.h"
#include <assert.h>
#include <math.h>

void insert_random_loops(struct tangle_state *tangle, int N)
{
  //TODO: configurable range of loop radii?, number of points?
  double D = 1;
  double Ls[3];
  for(int k = 0; k<3; ++k) {
	  Ls[k] = tangle->box.top_right_front.p[k] - tangle->box.bottom_left_back.p[k];
	  if(Ls[k] < D)
	    D = Ls[k];
  }

  const double rmin = 0.05*D;
  const double rmax = 0.25*D;

  for(int k = 0; k<N; ++k) {
    struct vec3d dir = vec3(drand48()-0.5, drand48()-0.5, drand48()-0.5);
    vec3_normalize(&dir);

    struct vec3d c;
    for(int j=0; j<3; ++j)
	    c.p[j] = tangle->box.bottom_left_back.p[j] + drand48()*Ls[j];

    double r =  rmin + drand48()*(rmax - rmin);
    add_circle(tangle, &c, &dir, r, 256);
  }
}

void insert_wall_bound_loops(struct tangle_state *tangle, int N)
{
  //assuming that the walls are in xy planes at zmin and zmax
  double D = 1;
  double Ls[3];
  for(int k = 0; k<3; ++k)
	  Ls[k] = tangle->box.top_right_front.p[k] - tangle->box.bottom_left_back.p[k];
  
  D = Ls[2];

  const double r = 0.5*D;

  for(int k = 0; k<N; ++k) {
    struct vec3d dir = vec3(drand48()-0.5, drand48()-0.5, drand48()-0.5);
    vec3_normalize(&dir);

    struct vec3d c;
    for(int j=0; j<2; ++j)
	    c.p[j] = tangle->box.bottom_left_back.p[j] + drand48()*Ls[j];

    if(drand48() < 0.5)
      c.p[2] = tangle->box.bottom_left_back.p[2];
    else
      c.p[2] = tangle->box.top_right_front.p[2];

    add_wall_circle(tangle, &c, &dir, r, 256);
  }
}

void make_big_ring(struct tangle_state *tangle, double ring_r, int ring_N)
{
  /*
   * Creates a big vortex ring at [0,0,0] with thickness of 10% of radius
   */

  struct vec3d dir = vec3(0, 0, 1);


  for(int k=0; k<ring_N; ++k) {
    //get some randomness in [-1, 1]
    double sx = (drand48() - 0.5)*2;
    double sy = (drand48() - 0.5)*2;
    double sz = (drand48() - 0.5)*2;

    struct vec3d cent = vec3(sx,sy,sz);
    vec3_mul(&cent, &cent, 0.1*ring_r);

    double sr = (drand48() - 0.5)*2;
    double r = (1 + sr/10)*ring_r;

    add_circle(tangle, &cent, &dir, r, 64);
  }
}

void random_straight_lines(struct tangle_state *tangle, int npairs, int points_per_line)
{
  //for the slab geometry, populate the computational box with straight vortices
  assert((tangle->box.wall[Z_L] == WALL_MIRROR && tangle->box.wall[Z_H] == WALL_MIRROR) ||
         (tangle->box.wall[Z_L] == WALL_PERIODIC && tangle->box.wall[Z_H] == WALL_PERIODIC));

  const double xmin = tangle->box.bottom_left_back.p[0];
  const double xmax = tangle->box.top_right_front.p[0];
  const double ymin = tangle->box.bottom_left_back.p[1];
  const double ymax = tangle->box.top_right_front.p[1];

  for(int k=0; k < npairs; ++k) {
    double x1 = xmin + (xmax - xmin)*drand48();
    double y1 = ymin + (ymax - ymin)*drand48();
    double x2 = xmin + (xmax - xmin)*drand48();
    double y2 = ymin + (ymax - ymin)*drand48();
    add_line(tangle, x1, y1, +1, points_per_line);
    add_line(tangle, x2, y2, -1, points_per_line);
  }
}

void random_polarized_lines(struct tangle_state *tangle, int nvort, int direction, int points_per_line)
{
  assert((tangle->box.wall[Z_L] == WALL_MIRROR && tangle->box.wall[Z_H] == WALL_MIRROR) ||
         (tangle->box.wall[Z_L] == WALL_PERIODIC && tangle->box.wall[Z_H] == WALL_PERIODIC));

  const double xmin = tangle->box.bottom_left_back.p[0];
  const double xmax = tangle->box.top_right_front.p[0];
  const double ymin = tangle->box.bottom_left_back.p[1];
  const double ymax = tangle->box.top_right_front.p[1];

  for(int k=0; k < nvort; ++k) {
    double x1 = xmin + (xmax - xmin)*drand48();
    double y1 = ymin + (ymax - ymin)*drand48();
    add_line(tangle, x1, y1, direction, points_per_line);
  }
}

void triangular_lattice(struct tangle_state *tangle, int shells, int direction, int points_per_line, float lattice_constant)
{
  const double dy = sqrt(3)/2;

  for(int i=-shells; i<shells+1; ++i)
  {
    for(int j=-shells; j<shells+1; ++j)
    {
      if(abs(i-j) > shells)
        continue;
      double x = i - 0.5*j;
      double y = dy*j;
      add_line(tangle, x*lattice_constant, y*lattice_constant, direction, points_per_line);
    }
  }
}

void quad_straight_lines(struct  tangle_state *tangle, int points_per_line, double separation)
{
  assert(tangle->box.wall[Z_L] == WALL_MIRROR && tangle->box.wall[Z_H] == WALL_MIRROR);

  add_line(tangle, separation/2, separation/2, +1, points_per_line);
  add_line(tangle, -separation/2, -separation/2, +1, points_per_line);
  add_line(tangle, -separation/2, separation/2, -1, points_per_line);
  add_line(tangle, separation/2, -separation/2, -1, points_per_line);
}

void dipole_straight_lines(struct  tangle_state *tangle, int points_per_line, int direction, double separation)
{
  assert(tangle->box.wall[Z_L] == WALL_MIRROR && tangle->box.wall[Z_H] == WALL_MIRROR);
 
  if(direction == 0) {
    add_line(tangle, separation/2, 0, +1, points_per_line);
    add_line(tangle, -separation/2, 0, -1, points_per_line);
  } else {
    add_line(tangle, 0, separation/2, +1, points_per_line);
    add_line(tangle, 0, -separation/2, -1, points_per_line);
  }
}

void single_straight_line(struct tangle_state *tangle, int points_per_line, int direction)
{
  assert(tangle->box.wall[Z_L] == WALL_MIRROR && tangle->box.wall[Z_H] == WALL_MIRROR);
  add_line(tangle, 0, 0, direction, points_per_line);
}