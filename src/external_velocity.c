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

#include <external_velocity.h>
#include <math.h>
#include <string.h>
#include "vec3_maths.h"
#include <stdio.h>

/*
 * Calculates the normal fluid velocity at point *where and
 * saves it into *res.
 *
 * returns 1 on success
 */

int get_vn(const struct vec3d *where, double t, struct vec3d *res)
{
  vn_conf.fun(where, t, res, &vn_conf);
  return 1;
}

int get_vs(const struct vec3d *where, double t, struct vec3d *res)
{
  vs_conf.fun(where, t, res, &vs_conf);
  return 1;
}

int get_vb(const struct vec3d *where, double t, struct vec3d *res)
{
  vb_conf.fun(where, t, res, &vb_conf);
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

struct v_conf_t vn_conf;
struct v_conf_t vs_conf;
struct v_conf_t vb_conf;

const struct v_conf_t v_conf_noflow = {
    .name = "no flow",
    .fun = get_v_noflow,
    .n_params = 0,
    .v_params = {{0,0,{.scalar=0}}}
};

struct v_conf_t v_confs[] = {
  {
    .name = "no flow",
    .fun = get_v_noflow,
    .n_params = 0,
    .v_params = {{0,0,{.scalar = 0}}}
  },
  {
    .name = "spherical",
    .fun = get_v_spherical,
    .n_params = 2,
    .v_params = {
      {"strength", scalar_param, {.scalar = 0}},
      {"cutoff", scalar_param, {.scalar = 0}}
    }
  },
  {
    .name = "simple",
    .fun = get_v_simple,
    .n_params = 1,
    .v_params = {
      {"external_v", vector_param, {.scalar = 0}}
    }
  },
  {
    .name = "spherical and simple",
    .fun = get_v_spherical_and_simple,
    .n_params = 3,
    .v_params = {
      {"external_v", vector_param, {.scalar = 0}},
      {"strength", scalar_param, {.scalar = 0}},
      {"cutoff", scalar_param, {.scalar = 0}}
    }
  },
  {
	.name = "cylindrical",
	.fun = get_v_cylindrical,
	.n_params = 2,
	.v_params = {
	    {"strength", scalar_param, {.scalar = 0}},
	    {"cutoff", scalar_param, {.scalar = 0}}
	}
  },
  {
    .name = "oscillating",
    .fun = get_v_oscillating,
    .n_params = 4,
    .v_params = {
      {"strength", scalar_param, {.scalar = 0}},
      {"frequency", scalar_param, {.scalar = 0}},
      {"wave", vector_param, {.scalar = 0}},
      {"polarization", vector_param, {.scalar = 0}}
    }
  },
  {
    .name = "oscillating slit",
    .fun = get_v_oscillating_slit,
    .n_params = 3,
    .v_params = {
      {"strength", scalar_param, {.scalar = 0}},
      {"frequency", scalar_param, {.scalar = 0}},
      {"cutoff", scalar_param, {.scalar = 0}}
    }
  },
  {
    .name = "coscos",
    .fun = get_v_coscos,
    .n_params = 3,
    .v_params = {
      {"k", scalar_param, {.scalar = 0}},
      {"v0", scalar_param, {.scalar = 0}},
      {"l", scalar_param, {.scalar = 0}}
    }
  },
  {
    .name = "coscos-divfree",
    .fun = get_v_coscos_divfree,
    .n_params = 3,
    .v_params = {
      {"k", scalar_param, {.scalar = 0}},
      {"v0", scalar_param, {.scalar = 0}},
      {"l", scalar_param, {.scalar = 0}}
    }
  },
  {
    .name = "cos-divfree",
    .fun = get_v_cos_divfree,
    .n_params = 3,
    .v_params = {
      {"k", scalar_param, {.scalar = 0}},
      {"v0", scalar_param, {.scalar = 0}},
      {"l", scalar_param, {.scalar = 0}}
    }
  },
  {
    .name = "simple-shear",
    .fun = get_v_simple_shear,
    .n_params = 2,
    .v_params = {
      {"rate", scalar_param, {.scalar=0}},
      {"v0", scalar_param, {.scalar=0}}
    }
  },
  {
    .name = "rotosc",
    .fun = get_v_rotosc,
    .n_params = 2,
    .v_params = {
      {"freq", scalar_param, {.scalar=0}},
      {"amp", scalar_param, {.scalar=0}}
    }
  },
  {
    .name = "rotosc-boundary",
    .fun = get_v_rotosc,
    .n_params = 2,
    .v_params = {
      {"freq", scalar_param, {.scalar=0}},
      {"amp", scalar_param, {.scalar=0}},
      {"delta", scalar_param, {.scalar=0}}
    }
  },
  {
    .name = "constant acceleration",
    .fun = get_v_constant_acceleration,
    .n_params = 1,
    .v_params = {
      {"accel", scalar_param, {.scalar=0}}
    }
  },
  {
    .name = "",
    .fun = NULL,
    .n_params = 0,
    .v_params = {{0,0,{.scalar = 0}}}
  }
};

int get_v_noflow(const struct vec3d *where __attribute__((unused)), double t __attribute__((unused)), struct vec3d *res, struct v_conf_t * v_conf __attribute__((unused)))
{
  *res = vec3(0, 0, 0); //no flow

  return 0;
}

int get_v_simple(const struct vec3d *where __attribute__((unused)), double t __attribute__((unused)), struct vec3d *res,  struct v_conf_t *vconf)
{
  return get_v_param_vector(vconf, "external_v", res);
}

/*
 * Returns a spherical sink (negative strength) or source (positive strength) velocity field.
 * The singularity is at origin and the magnitude is strength/(4 pi r^2) where r is the distance
 * from the origin to he point of interest (*where).
 */

int get_v_spherical(const struct vec3d *where, double t __attribute__((unused)), struct vec3d *res, struct v_conf_t *vconf)
{
  double strength;
  double cutoff;
  int err;
  if(!(err = get_v_param_scalar(vconf, "strength", &strength)))
    return err;
  if(!(err = get_v_param_scalar(vconf, "cutoff", &cutoff)))
    return err;

  double r = vec3_d(where);
  double attn = (cutoff/r)*(cutoff/r);
  if(r < 1e-8 || attn > 100) //hard cutoff to prevent 1/0 or underflow SIGFPE, probably not portable
    {
      *res = vec3(0,0,0);
      return 0;
    }

  double factor = exp(-attn);

  *res = *where;
  vec3_normalize(res);
  vec3_mul(res, res, strength/(4*M_PI*r*r)*factor);

  return 0;
}

int get_v_spherical_and_simple(const struct vec3d *where, double t __attribute__((unused)), struct vec3d *res, struct v_conf_t *vconf)
{
  struct vec3d v_spherical, v_simple;

  int err;
  err = get_v_spherical(where, t, &v_spherical, vconf);
  if(err) return err;

  err = get_v_simple(where, t, &v_simple, vconf);
  if(err) return err;

  vec3_add(res, &v_spherical, &v_simple);
  return 0;
}

int get_v_cylindrical(const struct vec3d *where, double t __attribute__((unused)), struct vec3d *res, struct v_conf_t *vconf)
{
  //strength is meant per unit length
  double strength;
  double cutoff;
  int err;
  if(!(err = get_v_param_scalar(vconf, "strength", &strength)))
    return err;
  if(!(err = get_v_param_scalar(vconf, "cutoff", &cutoff)))
    return err;

  //the cylinder is oriented along the z-axis
  double r = sqrt(where->p[0]*where->p[0] + where->p[1]*where->p[1]);
  double attn = (cutoff/r)*(cutoff/r);

  if(r < 1e-8 || attn > 100)
    {
      *res = vec3(0, 0, 0);
      return 0;
    }
  double factor = exp(-attn);

  *res = *where;
  res->p[2] = 0;
  vec3_normalize(res);
  vec3_mul(res, res, strength/(2*M_PI*r)*factor);
  return 0;
}

int get_v_oscillating(const struct vec3d *where, double t, struct vec3d *res, struct v_conf_t *vconf)
{
  double freq;
  double strength;
  struct vec3d k;
  struct vec3d pol;

  int err;
  if(!(err = get_v_param_scalar(vconf, "frequency", &freq)))
    return err;
  if(!(err = get_v_param_scalar(vconf, "strength", &strength)))
    return err;
  if(!(err = get_v_param_vector(vconf, "wave", &k)))
    return err;
  if(!(err = get_v_param_vector(vconf, "polarization", &pol)))
    return err;

  double X = cos(2*M_PI*vec3_dot(&k, where));
  double T = sin(2*M_PI*freq*t);

  vec3_mul(res, &pol, strength*X*T);
  return 0;
}

int get_v_oscillating_slit(const struct vec3d *where, double t, struct vec3d *res, struct v_conf_t *vconf)
{
  double freq;
  double strength;
  double cutoff;

  int err;
  if(!(err = get_v_param_scalar(vconf, "frequency", &freq)))
    return err;
  if(!(err = get_v_param_scalar(vconf, "strength", &strength)))
    return err;
  if(!(err = get_v_param_scalar(vconf, "cutoff", &cutoff)))
      return err;

  //slit is along y
  double x = where->p[0];
  double z = where->p[2];
  double r = sqrt(x*x + z*z);
  struct vec3d dir = vec3(x/r, 0, z/r);

  double attn = exp(-cutoff*cutoff/r/r);
  double osc = cos(2*M_PI*freq*t);

  vec3_mul(res, &dir, osc*attn);

  return 0;
}

int get_v_coscos(const struct vec3d *where, double t __attribute__((unused)), struct vec3d *res, struct v_conf_t *vconf)
{
  double k, v0, l;

  double x = where->p[0];
  double y = where->p[1];
  double z = where->p[2];

  int err;
  if(!(err = get_v_param_scalar(vconf, "k", &k)))
    return err;
  if(!(err = get_v_param_scalar(vconf, "v0", &v0)))
    return err;
  if(!(err = get_v_param_scalar(vconf, "l", &l)))
    return err;

  *res = vec3(0,0,0);
  res->p[2] = v0*(1 + cos(k*x)*cos(k*y)*exp(-z/l));

  return 0;
}

/*divergence-free coscos flow field*/
int get_v_coscos_divfree(const struct vec3d *where, double t __attribute__((unused)), struct vec3d *res, struct v_conf_t *vconf)
{
  double k, v0, l;

  double x = where->p[0];
  double y = where->p[1];
  double z = where->p[2];

  int err;
  if(!(err = get_v_param_scalar(vconf, "k", &k)))
    return err;
  if(!(err = get_v_param_scalar(vconf, "v0", &v0)))
    return err;
  if(!(err = get_v_param_scalar(vconf, "l", &l)))
    return err;

  *res = vec3(0,0,0);
  res->p[0] = v0/2/l/k*sin(k*x)*cos(k*y)*exp(-z/l);
  res->p[1] = v0/2/l/k*cos(k*x)*sin(k*y)*exp(-z/l);
  res->p[2] = v0*(1 + cos(k*x)*cos(k*y)*exp(-z/l));

  return 0;
}

int get_v_cos_divfree(const struct vec3d *where, double t __attribute__((unused)), struct vec3d *res, struct v_conf_t *vconf)
{
  double k, v0, l;

  double x = where->p[0];
  double z = where->p[2];

  int err;
  if(!(err = get_v_param_scalar(vconf, "k", &k)))
    return err;
  if(!(err = get_v_param_scalar(vconf, "v0", &v0)))
    return err;
  if(!(err = get_v_param_scalar(vconf, "l", &l)))
    return err;

  *res = vec3(0,0,0);
  res->p[0] = v0/l/k*sin(k*x)*exp(-z/l);
  res->p[1] = 0;
  res->p[2] = v0*(1 + cos(k*x)*exp(-z/l));

  return 0;
}

int get_v_simple_shear(const struct vec3d *where, double t __attribute__((unused)), struct vec3d *res, struct v_conf_t *vconf)
{
  double rate;
  double v0;
  double x = where->p[0];

  int err;
  if(!(err = get_v_param_scalar(vconf, "rate", &rate)))
    return err;
  if(!(err = get_v_param_scalar(vconf, "v0", &v0)))
      return err;

  *res = vec3(0,0,0);
  res->p[2] = v0 + rate*x;

  return 0;
}

int get_v_rotosc(const struct vec3d *where, double t, struct vec3d *res, struct v_conf_t *vconf)
{
  double freq;
  double amp;

  int err;
  if(!(err = get_v_param_scalar(vconf, "freq", &freq)))
    return err;
  if(!(err = get_v_param_scalar(vconf, "amp", &amp)))
    return err;

  const struct vec3d axis = vec3(0, 0, 1);
  struct vec3d v;
  vec3_cross(&v, &axis, where);
  vec3_mul(&v, &v, amp*cos(2*M_PI*freq*t));

  *res = v;
  return 0;
}

int get_v_rotosc_boundary(const struct vec3d *where, double t, struct vec3d *res, struct v_conf_t *vconf)
{
  double freq;
  double amp;
  double delta;

  int err;
  if(!(err = get_v_param_scalar(vconf, "freq", &freq)))
    return err;
  if(!(err = get_v_param_scalar(vconf, "amp", &amp)))
    return err;
  if(!(err = get_v_param_scalar(vconf, "delta", &delta)))
      return err;

  const struct vec3d axis = vec3(0, 0, 1);
  struct vec3d v;
  vec3_cross(&v, &axis, where);
  double z = where->p[2];
  vec3_mul(&v, &v, amp*exp(-z/delta)*cos(freq*t - z/delta));

  *res = v;
  return 0;
}

int get_v_constant_acceleration(const struct vec3d *where __attribute__((unused)), double t, struct vec3d *res, struct v_conf_t *vconf)
{
  double acc;
  int err;
  if(!(err = get_v_param_scalar(vconf, "accel", &acc)))
    return err;
  
  struct vec3d v = vec3(acc*t, 0, 0);
  *res = v;
  return 0;
}