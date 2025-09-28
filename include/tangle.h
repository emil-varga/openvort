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
extern const struct boundary_images periodic_z_open_xy;

extern const struct boundary_images wall_1_open;
extern const struct boundary_images wall_1_6;
extern const struct boundary_images wall_1_18;
extern const struct boundary_images wall_1_26;

extern const struct boundary_images wall_2_4;
extern const struct boundary_images wall_2_2;
extern const struct boundary_images wall_2_26;

extern const struct boundary_images channel_z;

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
  double *dxi; //arc length segments

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

struct octree;
/// @brief Update velocities at all nodes in the tangle.
/// @param tangle The vortex tangle, tangents and normals must be valid before calling this
/// @param t time for external velocities and boundary conditions, if used
/// @param tree octree used in the Barnes-Hut approximation, used only if the global use_BH is non-zero
///             if NULL and use_BH is True, new tree is built internally and destroyed at the end 
void update_velocities(struct tangle_state *tangle, double t, struct octree *tree);


/// @brief Update velocity at the node k of the tangle, this is called by update_valocities
/// @param tangle same as update_velocities
/// @param k index of the node to update
/// @param t same as update_velocities
/// @param tree If non-NULL and global use_BH is true, use the Barnes-Hut approximation represented by tree
void update_velocity(struct tangle_state *tangle, int k, double t, struct octree *tree);


/// @brief calculates the superfluid velocity
/// Calculates and returns the superfluid velocity v_s induced by tangle at location r.
/// Optionally disables one node given by skip. Use skip=-1 to not skip anything.
struct vec3d calculate_vs(struct tangle_state *tangle, struct vec3d r, int skip);

/// @brief Integration of the Biot-Savart integral
/// @param tangle The vortex tangle
/// @param r The point at which we want the velocity
/// @param skip Which index to skip (-1 to disable)
/// @param shift Calculate the velocity at r + shift (NULL to disable)
/// @param use_only_points Only use these points in the tangle (NULL to disable)
/// @param Npoints Number of points in use_only_points (ignored if use_only_points is NULL)
/// @return velocity at point *r*
struct vec3d calculate_vs_shift(const struct tangle_state *tangle, struct vec3d r, int skip, const struct vec3d *shift,
                                const int *use_only_points, const int Npoints);

void remesh(struct tangle_state *tangle, double min_dist, double max_dist);
void eliminate_small_loops(struct tangle_state *tangle, int loop_length);

int get_tangle_next_free(struct tangle_state *tangle);

static inline void update_tangle(struct tangle_state *tangle, double t)
{
  update_tangents_normals(tangle);
  update_velocities(tangle, t, NULL);
}

void remove_point(struct tangle_state *tangle, int point_idx, int merge_direction);
int add_point(struct tangle_state *tangle, int p);

int curvature_smoothing(struct tangle_state *tangle, double max_spp, double damping);
#endif //TANGLE_H
