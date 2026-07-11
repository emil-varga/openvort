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

#ifndef VEC3_MATHS_H
#define VEC3_MATHS_H

struct vec3d {
  double p[3];
};

struct mat3d {
  double m[3][3];
};

extern const double epsilon[3][3][3];
extern const struct mat3d mat_identity;

struct segment {
  struct vec3d r1, r2;
};

typedef enum _wt {WALL_OPEN, WALL_PERIODIC, WALL_MIRROR} wall_type;

/*
 * Faces of the computational box. x low, x high etc.
 *
 * Functions dealing with boundaries in tangle.c depend on
 * the ordering of this enum.
 */
typedef enum _fc {
  NOT_A_FACE = -1,
  X_L, X_H,
  Y_L, Y_H,
  Z_L, Z_H
} boundary_faces;

/*
 * Inward-facing normals of the box boundary face walls.
 * Can (and should) be indexed with boundary_faces enum.
 *
 * Defined in vec3_maths.c
 */
extern const struct vec3d boundary_normals[6];

struct domain_box {
  // Bottom/top refers to z (increasing from top to bottom)
  // Left/right refers to y (increasing from left to right)
  // Back/front refers to x (increasing from back to front)
  struct vec3d bottom_left_back;
  struct vec3d top_right_front;
  wall_type wall[6];
};

extern const struct vec3d DIR_X, DIR_Y, DIR_Z;
extern const struct vec3d DIRS[3];

#include <math.h>

static inline struct vec3d vec3(double x, double y, double z) {
  struct vec3d out = {.p = {x, y, z}};
  return out;
}

static inline void vec3_assign(struct vec3d *v, double x, double y, double z) {
  v->p[0] = x;
  v->p[1] = y;
  v->p[2] = z;
}

struct domain_box make_box(struct vec3d bottom_left_front,
			   struct vec3d top_right_back,
			   wall_type wall[6]);

int in_box(const struct domain_box *box, const struct vec3d *vec);
double max_box_size(const struct domain_box *box);
double min_box_size(const struct domain_box *box);
/*
 * open-space geometry
 */

static inline double vec3_dot(const struct vec3d *u, const struct vec3d *v) {
  double out = 0;
  for (int k = 0; k < 3; ++k)
    out += u->p[k] * v->p[k];

  return out;
}

static inline double vec3_d(const struct vec3d *u) {
  double out = 0;

  for (int k = 0; k < 3; ++k)
    out += u->p[k] * u->p[k];

  return sqrt(out);
}

//normalized dot, cosine of the angle
static inline double vec3_ndot(const struct vec3d *u, const struct vec3d *v) {
  double dot = vec3_dot(u, v);
  double d1 = vec3_d(u);
  double d2 = vec3_d(v);

  return dot / d1 / d2;
}

static inline void vec3_cross(struct vec3d *res,
		const struct vec3d *u, const struct vec3d *v) {
  struct vec3d vv = *v;
  struct vec3d uu = *u;
  res->p[0] = uu.p[1] * vv.p[2] - uu.p[2] * vv.p[1];
  res->p[1] = uu.p[2] * vv.p[0] - uu.p[0] * vv.p[2];
  res->p[2] = uu.p[0] * vv.p[1] - uu.p[1] * vv.p[0];
}

static inline void vec3_mul(struct vec3d *res,
	      const struct vec3d *u, double m) {
  for (int k = 0; k < 3; ++k)
    res->p[k] = u->p[k] * m;
}

static inline void vec3_sub(struct vec3d *res,
	      const struct vec3d *u, const struct vec3d *v) {
  for (int k = 0; k < 3; ++k)
    res->p[k] = u->p[k] - v->p[k];
}

static inline void vec3_add(struct vec3d *res,
	      const struct vec3d *u, const struct vec3d *v) {
  for (int k = 0; k < 3; ++k)
    res->p[k] = u->p[k] + v->p[k];
}

static inline struct vec3d vec3_add2(const struct vec3d *u, const struct vec3d *v) {
  struct vec3d out;
  vec3_add(&out, u, v);

  return out;
}

static inline double vec3_dist(const struct vec3d *u, const struct vec3d *v) {
  struct vec3d vv;
  vec3_sub(&vv, u, v);
  return vec3_d(&vv);
}

static inline void vec3_normalize(struct vec3d *v) {
  double d = vec3_d(v);

  vec3_mul(v, v, 1 / d);
}

static inline struct vec3d segment_to_vec(const struct segment *seg) {
  struct vec3d v;
  vec3_sub(&v, &seg->r2, &seg->r1);

  return v;
}

static inline double segment_len(const struct segment *seg)
{
  return vec3_dist(&seg->r1, &seg->r2);
}

/*
 * Periodic and mirror geometries
 */

struct segment seg_pwrap(const struct vec3d *r1, const struct vec3d *r2,
			 const struct domain_box *box);

struct vec3d box_shift(const struct vec3d *v, const struct domain_box *box,
		       int shift[3]);

struct vec3d mirror_shift(const struct vec3d *v, const struct domain_box *box,
			  boundary_faces wall);
struct vec3d mirror_dir_reflect(const struct vec3d *v, boundary_faces wall);

struct vec3d periodic_shift(const struct vec3d *v, const struct domain_box *box,
			    boundary_faces wall);

/*
 * Matrix operations
 */

void mat3_dmul(struct mat3d *res, const struct mat3d *M, double m);
struct mat3d mat3_null();
void mat3_add(struct mat3d *res, const struct mat3d *a, const struct mat3d *b);
void mat3_sub(struct mat3d *res, const struct mat3d *a, const struct mat3d *b);
void mat3_mul(struct mat3d *res, const struct mat3d *a, const struct mat3d *b);
void mat3_vmul(struct vec3d *res, const struct mat3d *a, const struct vec3d *v);
void vec3_outer(struct mat3d *res, const struct vec3d *u, const struct vec3d *v);//outer product

#endif//VEC3_MATHS_H
