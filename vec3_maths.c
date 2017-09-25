#include "vec3_maths.h"
#include <math.h>
#include <stdio.h>

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
  res->p[0] = u->p[1] * v->p[2] - u->p[2] * v->p[1];
  res->p[1] = u->p[2] * v->p[0] - u->p[0] * v->p[2];
  res->p[2] = u->p[0] * v->p[1] - u->p[1] * v->p[0];
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

double vec3_d(const struct vec3d *u)
{
  double out = 0;

  for(int k=0; k<3; ++k)
    out += u->p[k]*u->p[k];

  return sqrt(out);
}

double vec3_dist(const struct vec3d *u, const struct vec3d *v)
{
  struct vec3d vv;
  vec3_sub(&vv, u, v);
  return vec3_d(&vv);
}
