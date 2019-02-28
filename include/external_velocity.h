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

#ifndef EXTERNAL_VELOCITY_H
#define EXTERNAL_VELOCITY_H

/* Here we deal with setting up and calculating externally-imposed velocity fields.
 * Both non-homogeneneous and non-stationary flows are possible.
 *
 * The flow is selected through the configuration file with the parameters "vn_conf" and
 * "vs_conf" (see configuration.c). The parameters of the flow are declared through the
 * v_confs[] array and specified through the configuration file.
 *
 * The selected flow field is available through get_vn and get_vs functions.
 *
 * New flows should be added to external_velocity.c and the flow function should be
 * declared here (see below).
 */

#include "vec3_maths.h"

/*
 * These return the external normal or superfluid velocity at point *where at time t.
 * The velocity is retured in *res. Return 0 on success.
 */

int get_vn(const struct vec3d *where, double t, struct vec3d *res);
int get_vs(const struct vec3d *where, double t, struct vec3d *res);
int get_vb(const struct vec3d *where, double t, struct vec3d *res);

/*
 * struct v_param_t holds information about a single parameter of the flow.
 *
 * The parameter has a name which is the same through which it is set in the configuration file.
 * Parameters can be either scalars or 3-vectors.
 */

#define scalar_param 1
#define vector_param 2

struct v_param_t {
  char *name;
  int type;
  union param_value_t {
    struct vec3d vector;
    double scalar;
  } value;
};


/*
 * Function prototype for the flow field function. The function should take position,
 * time, return the result in the 3rd parameter and should look for its parameters in the
 * vconf array which will be supplied to it.
 */

struct v_conf_t;
typedef int (*get_v_t)(const struct vec3d *, double, struct vec3d *, struct v_conf_t *vconf);

/*
 * Configuration of the flow.
 */

struct v_conf_t {
  char *name;    //the name of the configuration
  get_v_t fun;   //function which will do the actual calculation
  int n_params;  //how many parameters does the flow conf have (at most 32)
  struct v_param_t v_params[32]; //the array of parameters, at most 32.
};

extern struct v_conf_t v_confs[];   //all flow confs
extern struct v_conf_t vn_conf;     //normal fluid conf, will be taken from v_confs[] during configuration
extern struct v_conf_t vs_conf;     //superfluid conf, will be taken from v_confs[] during configuration
extern struct v_conf_t vb_conf; //the moving boundary only works for the bottom z-wall

//helper functions to look up scalar or vector parameters in the v_conf_t* parameter arrays
int get_v_param_scalar(const struct v_conf_t *vconf, const char *name, double *res);
int get_v_param_vector(const struct v_conf_t *vconf, const char *name, struct vec3d *res);

/*
 * Currently available flow configurations.
 */
int get_v_noflow(const struct vec3d *where, double t, struct vec3d *res, struct v_conf_t *vconf);
int get_v_simple(const struct vec3d *where, double t, struct vec3d *res, struct v_conf_t *vconf);
int get_v_spherical(const struct vec3d *where, double t, struct vec3d *res, struct v_conf_t *vconf);
int get_v_spherical_and_simple(const struct vec3d *where, double t, struct vec3d *res, struct v_conf_t *vconf);
int get_v_cylindrical(const struct vec3d *where, double t, struct vec3d *res, struct v_conf_t *vconf);
int get_v_oscillating(const struct vec3d *where, double t, struct vec3d *res, struct v_conf_t *vconf);
int get_v_coscos(const struct vec3d *where, double t, struct vec3d *res, struct v_conf_t *vconf);
int get_v_coscos_divfree(const struct vec3d *where, double t, struct vec3d *res, struct v_conf_t *vconf);
int get_v_cos_divfree(const struct vec3d *where, double t, struct vec3d *res, struct v_conf_t *vconf);
int get_v_simple_shear(const struct vec3d *where, double t, struct vec3d *res, struct v_conf_t *vconf);

//oscillatory rotation, with attenuation and the boundary version with viscous penetratio depth
//the boundary version requires that the boundary is at 0
int get_v_rotosc(const struct vec3d *where, double t, struct vec3d *res, struct v_conf_t *vconf);
int get_v_rotosc_boundary(const struct vec3d *where, double t, struct vec3d *res, struct v_conf_t *vconf);

#endif /* EXTERNAL_VELOCITY_H */
