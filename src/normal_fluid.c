/*
 * normal_fluid.c
 *
 *  Created on: Jan 10, 2018
 *      Author: emil
 */

#include <math.h>
#include "normal_fluid.h"
#include "vec3_maths.h"

/*
 * Calculates the normal fluid velocity at point *where and
 * saves it into *res.
 *
 * returns 1 on success
 */
get_vn_t get_vn;

struct vn_conf vn_confs[] = {
    {
	.name = "no flow",
	.fun = get_vn_noflow,
	.n_params = 0,
	.vn_params = {{0,0,0}}
    },
    {
	.name = "spherical",
	.fun = get_vn_spherical,
	.n_params = 2,
	.vn_params = {
	    {"spherical_vn_strength", &spherical_vn_strength, scalar_param},
	    {"spherical_vn_cutoff", &spherical_vn_cutoff, scalar_param}
	}
    },
    {
	.name = "simple",
	.fun = get_vn_simple_cf,
	.n_params = 1,
	.vn_params = {
	    {"external_vn", &external_vn, vector_param}
	}
    }
};

int get_vn_noflow(const struct vec3d *where, struct vec3d *res)
{
  *res = vec3(0, 0, 0); //no flow

  return 0;
}

struct vec3d external_vn = {{0, 0, 0}};

void set_external_vn(double vx, double vy, double vz)
{
  external_vn = vec3(vx, vy, vz);
}

int get_vn_simple_cf(const struct vec3d *where, struct vec3d *res)
{
  *res = external_vn;

  return 0;
}

/*
 * Returns a spherical sink (negative strength) or source (positive strength) velocity field.
 * The singularity is at origin and the magnitude is strength/(4 pi r^2) where r is the distance
 * from the origin to he point of interest (*where).
 */

double spherical_vn_strength;
double spherical_vn_cutoff;
int get_vn_spherical(const struct vec3d *where, struct vec3d *res)
{
  double strength = spherical_vn_strength;
  double cutoff = spherical_vn_cutoff;
  double r = vec3_d(where);
  if(r < cutoff)
    {
      *res = vec3(0,0,0);
      return 0;
    }

  *res = *where;
  vec3_normalize(res);
  vec3_mul(res, res, strength/(4*M_PI*r*r));

  return 0;
}
