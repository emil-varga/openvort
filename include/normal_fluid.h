/*
 * normal_fluid.h
 *
 *  Created on: Jan 10, 2018
 *      Author: emil
 */

#ifndef INCLUDE_NORMAL_FLUID_H_
#define INCLUDE_NORMAL_FLUID_H_

#include "vec3_maths.h"

int get_vn(struct vec3d *where, struct vec3d *res);

extern struct vec3d external_vn;

#endif /* INCLUDE_NORMAL_FLUID_H_ */
