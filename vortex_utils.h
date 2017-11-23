#ifndef VORTEX_UTILS_H
#define VORTEX_UTILS_H

//This file contains miscellaneous utilities related to
//the handling of the tangle. Actual calculationg belong to tangle.h.

//This contains routines for creating basic vortex geometries and
//handling file I/O

#include <stdlib.h>
#include <stdio.h>
#include "vec3_maths.h"
#include "tangle.h"

void add_circle(struct tangle_state *tangle,
		struct vec3d *center, struct vec3d *dir, double r,
		size_t Npoints);

void save_tangle(const char *filename, struct tangle_state *tangle);

int check_integrity(const struct tangle_state *tangle);

#endif
