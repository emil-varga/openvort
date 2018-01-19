/*
 * normal_fluid.c
 *
 *  Created on: Jan 10, 2018
 *      Author: emil
 */

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
