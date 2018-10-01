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

#include "vec3_maths.h"

int get_vn(const struct vec3d *where, double t, struct vec3d *res);
int get_vs(const struct vec3d *where, double t, struct vec3d *res);

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

struct v_conf_t;
typedef int (*get_v_t)(const struct vec3d *, double, struct vec3d *, struct v_conf_t *vconf);

struct v_conf_t {
  char *name;
  get_v_t fun;
  int n_params;
  struct v_param_t v_params[32];
};

extern struct v_conf_t v_confs[];
extern struct v_conf_t vn_conf;
extern struct v_conf_t vs_conf;

int get_v_param_scalar(const struct v_conf_t *vconf, const char *name, double *res);
int get_v_param_vector(const struct v_conf_t *vconf, const char *name, struct vec3d *res);

int get_v_noflow(const struct vec3d *where, double t, struct vec3d *res, struct v_conf_t *vconf);
int get_v_simple(const struct vec3d *where, double t, struct vec3d *res, struct v_conf_t *vconf);
int get_v_spherical(const struct vec3d *where, double t, struct vec3d *res, struct v_conf_t *vconf);
int get_v_spherical_and_simple(const struct vec3d *where, double t, struct vec3d *res, struct v_conf_t *vconf);
int get_v_cylindrical(const struct vec3d *where, double t, struct vec3d *res, struct v_conf_t *vconf);
int get_v_oscillating(const struct vec3d *where, double t, struct vec3d *res, struct v_conf_t *vconf);

#endif /* EXTERNAL_VELOCITY_H */
