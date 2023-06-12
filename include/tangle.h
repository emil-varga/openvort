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

#ifndef TANGLE_H
#define TANGLE_H

#include <stdlib.h>
#include <stdint.h>

#include "vec3_maths.h"

typedef enum _ns_e {
  EMPTY = -1, //no point here
  FREE, //ordinary vortex point
  PINNED, //pinned on the wall
  PINNED_SLIP //pinned, but can slip
} node_status_t;

//status of the node and on which wall it is pinned
typedef struct _ns {
  node_status_t status;
  boundary_faces pin_wall;
} node_status;

#define PERIODIC_X 1
#define PERIODIC_Y 1<<1
#define PERIDOIC_Z 1<<2

struct neighbour_t {
  int forward;
  int reverse;
};

/*
 * Structures describing the virtual tangles used for implementation of boundary conditions
 */

struct image_tangle {
  int shift[3]; //positive/negative number of shifts in the units of box size
  int reflect; //either a boundary_faces enum index or -1 for periodic
};

struct boundary_images {
  const struct image_tangle *images;
  int n;
};

//the actual image tangle configurations
//defined in boundary_images.c
extern const struct boundary_images open_boundaries;
extern const struct boundary_images periodic_6;
extern const struct boundary_images periodic_18;
extern const struct boundary_images periodic_26;

extern const struct boundary_images wall_1_open;
extern const struct boundary_images wall_1_6;
extern const struct boundary_images wall_1_18;
extern const struct boundary_images wall_1_26;

extern const struct boundary_images wall_2_4;
extern const struct boundary_images wall_2_2;

/*
 * The structure that holds all the tangle information
 */
struct tangle_state {
  struct vec3d *vnodes;
  struct vec3d *vnodes_new;
  struct vec3d *vels; //node velocities
  struct vec3d *vs; //superfluid velocity at the node
  struct vec3d *tangents;
  struct vec3d *normals;

  //flags that the properties of the points
  //need to be recalculated
  //currently used for not reconnecting twice in a single pass
  int *recalculate;

  struct neighbour_t *connections;

  node_status *status; //pinned, free, empty, etc.
  //the boundary conditions in the 6 cardinal directions
  struct domain_box box;
  struct boundary_images bimg;

  int N;
  int next_free;
  int total_free;
};

void create_tangle(struct tangle_state *tangle, size_t n);
void expand_tangle(struct tangle_state *tangle, size_t n);
void free_tangle(struct tangle_state *tangle);
struct vec3d step_node(const struct tangle_state *tangle, int i, int where);
int num_free_points(struct tangle_state *tangle);
int tangle_total_points(struct tangle_state *tangle);

//calculates the shifted r according to the image_tangle conf
struct vec3d shifted(const struct image_tangle *shift, const struct tangle_state *tangle,
		     const struct vec3d *r);
void enforce_boundaries(struct tangle_state *tangle);
void update_tangent_normal(struct tangle_state *tangle, size_t k);
void update_tangents_normals(struct tangle_state *tangle);
void update_velocities(struct tangle_state *tangle, double t);
void update_velocity(struct tangle_state *tangle, int k, double t);

/**
  @brief calculates the superfluid velocity

  Calculates and returns the superfluid velocity v_s induced by tangle at location r.
  Optionally disables one node given by skip. Use skip=-1 to not skip anything.
*/
struct vec3d calculate_vs(struct tangle_state *tangle, struct vec3d r, int skip);

void remesh(struct tangle_state *tangle, double min_dist, double max_dist);
void eliminate_small_loops(struct tangle_state *tangle, int loop_length);

int get_tangle_next_free(struct tangle_state *tangle);

static inline void update_tangle(struct tangle_state *tangle, double t)
{
  update_tangents_normals(tangle);
  update_velocities(tangle, t);
}

void remove_point(struct tangle_state *tangle, int point_idx);
int add_point(struct tangle_state *tangle, int p);

int curvature_smoothing(struct tangle_state *tangle, double max_spp, double damping);
#endif //TANGLE_H
