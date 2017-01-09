#ifndef VORTEX_UTILS_H
#define VORTEX_UTILS_H

#include <stdlib.h>
#include "vec3_maths.h"
#include "tangle.h"

void add_circle(struct tangle_state *tangle,
		struct vec3d *center, struct vec3d *dir, double r,
		size_t Npoints);

void save_tangle(FILE *stream, struct tangle_state *tangle);

#endif
