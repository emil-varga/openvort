/*
 * normal_fluid.h
 *
 *  Created on: Jan 10, 2018
 *      Author: emil
 */

#ifndef INCLUDE_EXTERNAL_VELOCITY_H_
#define INCLUDE_EXTERNAL_VELOCITY_H_

#include "vec3_maths.h"

int get_vn(const struct vec3d *where, struct vec3d *res);
int get_vs(const struct vec3d *where, struct vec3d *res);

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
typedef int (*get_v_t)(const struct vec3d *, struct vec3d *, struct v_conf_t *vconf);

struct v_conf_t {
  char *name;
  get_v_t fun;
  int n_params;
  struct v_param_t v_params[32];
};

extern struct v_conf_t v_confs[];
extern struct v_conf_t *vn_conf;
extern struct v_conf_t *vs_conf;

int get_v_param_scalar(const struct v_conf_t *vconf, const char *name, double *res);
int get_v_param_vector(const struct v_conf_t *vconf, const char *name, struct vec3d *res);

int get_v_noflow(const struct vec3d *where, struct vec3d *res, struct v_conf_t *vconf);
int get_v_simple(const struct vec3d *where, struct vec3d *res, struct v_conf_t *vconf);
int get_v_spherical(const struct vec3d *where, struct vec3d *res, struct v_conf_t *vconf);

#endif /* INCLUDE_EXTERNAL_VELOCITY_H_ */
