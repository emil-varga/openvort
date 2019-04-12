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

#include "vec3_maths.h"
#include <math.h>
#include <stdio.h>
#include "util.h"

const struct vec3d DIR_X = {{1, 0, 0}};
const struct vec3d DIR_Y = {{0, 1, 0}};
const struct vec3d DIR_Z = {{0, 0, 1}};
const struct vec3d DIRS[3] = {{{1, 0, 0}}, {{0, 1, 0}}, {{0, 0, 1}}};

/*
 * Inward-facing normals of the box boundary face walls
 */

const struct vec3d boundary_normals[] = {
    {{1, 0, 0}},  // X_L
    {{-1, 0, 0}}, // X_H
    {{0, 1, 0}},  // Y_L
    {{0, -1, 0}}, // Y_H
    {{0, 0, 1}},  // Z_L
    {{0, 0, -1}}  // Z_H
};

/*
 * Is the given vector vec inside the box?
 */
int in_box(const struct domain_box *box, const struct vec3d *vec)
{
  double vx = vec->p[0];
  double vy = vec->p[1];
  double vz = vec->p[2];

  int in_box;
  in_box =
      vx > box->bottom_left_back.p[0] && vx < box->top_right_front.p[0] &&
      vy > box->bottom_left_back.p[1] && vy < box->top_right_front.p[1] &&
      vz > box->bottom_left_back.p[2] && vz < box->top_right_front.p[2];
  return in_box;
}

struct segment seg_pwrap(const struct vec3d *r1, const struct vec3d *r2,
			 const struct domain_box *box)
{
  struct segment seg = {
      .r1 = *r1,
      .r2 = *r2
  };

  double Lx = box->top_right_front.p[0] - box->bottom_left_back.p[0];
  double Ly = box->top_right_front.p[1] - box->bottom_left_back.p[1];
  double Lz = box->top_right_front.p[2] - box->bottom_left_back.p[2];

  double dx = seg.r2.p[0] - seg.r1.p[0];
  double dy = seg.r2.p[1] - seg.r1.p[1];
  double dz = seg.r2.p[2] - seg.r1.p[2];

  if(box->wall[X_L] == WALL_PERIODIC)
    {
      if(dx > Lx/2)
	seg.r2.p[0] -= Lx;
      else if(dx < -Lx/2)
	seg.r2.p[0] += Lx;
    }
  if(box->wall[Y_L] == WALL_PERIODIC)
    {
      if(dy > Ly/2)
	seg.r2.p[1] -= Ly;
      else if(dy < -Ly/2)
	seg.r2.p[1] += Ly;
    }
  if(box->wall[Z_L] == WALL_PERIODIC)
    {
      if(dz > Lz/2)
	seg.r2.p[2] -= Lz;
      else if(dz < -Lz/2)
	seg.r2.p[2] += Lz;
    }

  return seg;
}

struct vec3d box_shift(const struct vec3d *v, const struct domain_box *box,
		       int shift[3])
{
  double Lx = box->top_right_front.p[0] - box->bottom_left_back.p[0];
  double Ly = box->top_right_front.p[1] - box->bottom_left_back.p[1];
  double Lz = box->top_right_front.p[2] - box->bottom_left_back.p[2];
  double Ls[] = {Lx, Ly, Lz};

  struct vec3d out = *v;

  for(int k = 0; k < 3; ++k)
    out.p[k] += Ls[k]*shift[k];

  return out;
}

struct vec3d periodic_shift(const struct vec3d *v, const struct domain_box *box,
			  boundary_faces wall)
{
  struct vec3d mv = *v;

  int coord;

  double Lx = box->top_right_front.p[0] - box->bottom_left_back.p[0];
  double Ly = box->top_right_front.p[1] - box->bottom_left_back.p[1];
  double Lz = box->top_right_front.p[2] - box->bottom_left_back.p[2];
  double Ls[] = {Lx, Ly, Lz};

  switch(wall)
  {
    case X_L: case X_H: coord=0; break;
    case Y_L: case Y_H: coord=1; break;
    case Z_L: case Z_H: coord=2; break;
    default:
      error("Bad wall %d", wall);
      coord=-1;
      break;
  }

  switch(wall)
  {
    case X_L:
    case Y_L:
    case Z_L:
      mv.p[coord] += Ls[coord];
      break;
    case X_H:
    case Y_H:
    case Z_H:
      mv.p[coord] -= Ls[coord];
      break;
    default:
      error("Bad wall %d", wall);
      break;
  }

  return mv;
}

struct vec3d mirror_shift(const struct vec3d *v, const struct domain_box *box,
			    boundary_faces wall)
{
  struct vec3d mv = *v;

  int coord;
  double wall_pos;

  switch(wall)
  {
    case X_L: case X_H: coord=0; break;
    case Y_L: case Y_H: coord=1; break;
    case Z_L: case Z_H: coord=2; break;
    default:
      error("Bad wall %d", wall);
      coord=-1;
      break;
  }

  switch(wall)
  {
    case X_L: wall_pos = box->bottom_left_back.p[0]; break;
    case X_H: wall_pos = box->top_right_front.p[0]; break;
    case Y_L: wall_pos = box->bottom_left_back.p[1]; break;
    case Y_H: wall_pos = box->top_right_front.p[1]; break;
    case Z_L: wall_pos = box->bottom_left_back.p[2]; break;
    case Z_H: wall_pos = box->top_right_front.p[2]; break;
    default:
      error("Bad wall %d", wall);
      wall_pos = -1;
      break;
  }

  mv.p[coord] = 2*wall_pos - mv.p[coord];

  return mv;
}

struct vec3d mirror_dir_reflect(const struct vec3d *v, boundary_faces wall)
{
  /*
   * Reflect a directional vector in a mirror (i.e., the component normal
   * to the mirror reverses).
   */
  struct vec3d mv = *v;

  double nc = vec3_dot(&mv, &boundary_normals[wall]);
  struct vec3d tmp;
  vec3_mul(&tmp, &boundary_normals[wall], -2*nc);
  vec3_add(&mv, &mv, &tmp);

  return mv;
}

struct vec3d vec3(double x, double y, double z)
{
  struct vec3d out = {.p = {x, y, z}};
  return out;
}

void vec3_assign(struct vec3d *v, double x, double y, double z)
{
  v->p[0] = x;
  v->p[1] = y;
  v->p[2] = z;
}

double vec3_dot(const struct vec3d *u, const struct vec3d *v)
{
  double out = 0;
  for(int k=0; k<3; ++k)
    out += u->p[k]*v->p[k];

  return out;
}

double vec3_ndot(const struct vec3d *u, const struct vec3d *v)
{
  double dot = vec3_dot(u, v);
  double d1 = vec3_d(u);
  double d2 = vec3_d(v);

  return dot/d1/d2;
}

void vec3_cross(struct vec3d *res,
		const struct vec3d *u, const struct vec3d *v)
{
  struct vec3d vv = *v;
  struct vec3d uu = *u;
  res->p[0] = uu.p[1] * vv.p[2] - uu.p[2] * vv.p[1];
  res->p[1] = uu.p[2] * vv.p[0] - uu.p[0] * vv.p[2];
  res->p[2] = uu.p[0] * vv.p[1] - uu.p[1] * vv.p[0];
}

void vec3_mul(struct vec3d *res,
	      const struct vec3d *u, double m)
{
  for(int k=0; k<3; ++k)
    res->p[k] = u->p[k]*m;
}

void vec3_sub(struct vec3d *res,
	      const struct vec3d *u, const struct vec3d *v)
{
  for(int k=0; k<3; ++k)
    res->p[k] = u->p[k] - v->p[k];
}

void vec3_add(struct vec3d *res,
	      const struct vec3d *u, const struct vec3d *v)
{
  for(int k=0; k<3; ++k)
    res->p[k] = u->p[k] + v->p[k];
}

struct vec3d vec3_add2(const struct vec3d *u, const struct vec3d *v)
{
  struct vec3d out;
  vec3_add(&out, u, v);

  return out;
}

double vec3_d(const struct vec3d *u)
{
  double out = 0;

  for(int k=0; k<3; ++k)
    {
      out += u->p[k]*u->p[k];
    }

  return sqrt(out);
}

double vec3_dist(const struct vec3d *u, const struct vec3d *v)
{
  struct vec3d vv;
  vec3_sub(&vv, u, v);
  return vec3_d(&vv);
}

void vec3_normalize(struct vec3d *v)
{
  double d = vec3_d(v);

  vec3_mul(v, v, 1/d);
}

struct vec3d segment_to_vec(const struct segment *seg)
{
  struct vec3d v;
  vec3_sub(&v, &seg->r2, &seg->r1);

  return v;
}

void mat3_dmul(struct mat3d *res, const struct mat3d *M, double m)
{
  for(int i=0; i<3; ++i)
    for(int j=0; j<3; ++j)
      res->m[i][j] = m*M->m[i][j];
}
struct mat3d mat3_null()
{
  struct mat3d m;
  for(int i=0; i<3; ++i)
    for(int j=0; j<3; ++j)
      m.m[i][j] = 0;
  return m;
}
void mat3_add(struct mat3d *res, const struct mat3d *a, const struct mat3d *b)
{
  for(int i=0; i<3; ++i)
    for(int j=0; j<3; ++j)
      res->m[i][j] = a->m[i][j] + b->m[i][j];
}
void mat3_sub(struct mat3d *res, const struct mat3d *a, const struct mat3d *b)
{
  for(int i=0; i<3; ++i)
    for(int j=0; j<3; ++j)
      res->m[i][j] = a->m[i][j] - b->m[i][j];
}
void mat3_mul(struct mat3d *res, const struct mat3d *_a, const struct mat3d *_b)
{
  //we need copies in case res points to the same as _a or _b
  struct mat3d a = *_a;
  struct mat3d b = *_b;
  for(int i=0; i<3; ++i)
    for(int j=0; j<3; ++j)
      for(int k=0; k<3; ++k)
	res->m[i][j] = a.m[i][k]*b.m[k][j];
}
void mat3_vmul(struct vec3d *res, const struct mat3d *a, const struct vec3d *v)
{
  for(int i=0; i<3; ++i)
    for(int j=0; j<3; ++j)
	res->p[i] = a->m[i][j]*v->p[j];
}
void vec3_outer(struct mat3d *res, const struct vec3d *u, const struct vec3d *v)
{
  //outer product
  for(int i=0; i<3; ++i)
      for(int j=0; j<3; ++j)
	  res->m[i][j] = u->p[i]*v->p[j];
}
