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
 * returns 0 on success
 *
 * TODO: now it just returns a constant, this should be generalized to more
 * 		 complicated flow
 */
int get_vn(struct vec3d *where, struct vec3d *res)
{
  *res = vec3(0, 0, 0); //no flow

  return 0;
}

struct vec3d external_vn = {0, 0, 0};

void set_external_vn(double vx, double vy, double vz)
{
  external_vn = vec3(vx, vy, vz);
}

int get_vn_simple_cf(struct vec3d *where, struct vec3d *res)
{
  *res = external_vn;

  return 0;
}

/*
 * Returns a spherical sink (negative strength) or source (positive strength) velocity field.
 * The singularity is at origin and the magnitude is strength/(4 pi r^2) where r is the distance
 * from the origin to he point of interest (*where).
 */
int get_vn_spherical(struct vec3d *where, struct vec3d *res, double strength, double cutoff)
{
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
