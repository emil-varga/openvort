/*
 * normal_fluid.h
 *
 *  Created on: Jan 10, 2018
 *      Author: emil
 */

#ifndef INCLUDE_NORMAL_FLUID_H_
#define INCLUDE_NORMAL_FLUID_H_

#include "vec3_maths.h"

typedef int (*get_vn_t)(const struct vec3d *, struct vec3d *);

extern get_vn_t get_vn;

#define scalar_param 1
#define vector_param 2

struct vn_param_t {
  const char *name;
  void *addr;
  const int type;
};

struct vn_conf {
  const char *name;
  get_vn_t fun;
  const int n_params;
  struct vn_param_t vn_params[32];
};

extern struct vn_conf vn_confs[];
int get_vn_noflow(const struct vec3d *where, struct vec3d *res);

extern struct vec3d external_vn;
int get_vn_simple_cf(const struct vec3d *where, struct vec3d *res);

extern double spherical_vn_strength;
extern double spherical_vn_cutoff;
int get_vn_spherical(const struct vec3d *where, struct vec3d *res);

#endif /* INCLUDE_NORMAL_FLUID_H_ */
