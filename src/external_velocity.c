/*
 * normal_fluid.c
 *
 *  Created on: Jan 10, 2018
 *      Author: emil
 */

#include <external_velocity.h>
#include <math.h>
#include <string.h>
#include "vec3_maths.h"

/*
 * Calculates the normal fluid velocity at point *where and
 * saves it into *res.
 *
 * returns 1 on success
 */

int get_vn(const struct vec3d *where, struct vec3d *res)
{
  vn_conf->fun(where, res, vn_conf);
  return 1;
}

int get_vs(const struct vec3d *where, struct vec3d *res)
{
  vs_conf->fun(where, res, vs_conf);
  return 1;
}

int get_v_param_scalar(const struct v_conf_t *vconf, const char *name, double *res)
{
  for(int k=0; k<vconf->n_params; ++k)
    {
      if(strcmp(name, vconf->v_params[k].name) == 0)
	{
	  *res = vconf->v_params[k].value.scalar;
	  return 1;
	}
    }

  return 0;
}

int get_v_param_vector(const struct v_conf_t *vconf, const char *name, struct vec3d *res)
{
  for(int k=0; k<vconf->n_params; ++k)
      {
        if(strcmp(name, vconf->v_params[k].name) == 0)
  	{
  	  *res = vconf->v_params[k].value.vector;
  	  return 1;
  	}
      }

  return 0;
}

struct v_conf_t *vn_conf = &v_confs[0];
struct v_conf_t *vs_conf = &v_confs[0];

struct v_conf_t v_confs[] = {
    {
        .name = "no flow",
        .fun = get_v_noflow,
        .n_params = 0,
        .v_params = {{0,0}}
    },
    {
	.name = "spherical",
	.fun = get_v_spherical,
	.n_params = 2,
	.v_params = {
	    {"strength", scalar_param},
	    {"cutoff", scalar_param}
	}
    },
    {
	.name = "simple",
	.fun = get_v_simple,
	.n_params = 1,
	.v_params = {
	    {"external_v", vector_param}
	}
    }
};

int get_v_noflow(const struct vec3d *where, struct vec3d *res, struct v_conf_t *vconf)
{
  *res = vec3(0, 0, 0); //no flow

  return 0;
}

int get_v_simple(const struct vec3d *where, struct vec3d *res,  struct v_conf_t *vconf)
{
  return get_v_param_vector(vconf, "simple", res);
}

/*
 * Returns a spherical sink (negative strength) or source (positive strength) velocity field.
 * The singularity is at origin and the magnitude is strength/(4 pi r^2) where r is the distance
 * from the origin to he point of interest (*where).
 */

int get_v_spherical(const struct vec3d *where, struct vec3d *res, struct v_conf_t *vconf)
{
  double strength;
  double cutoff;
  int err;
  if(err = get_v_param_scalar(vconf, "strength", &strength))
    return err;
  if(err = get_v_param_scalar(vconf, "cutoff", &cutoff))
    return err;

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
