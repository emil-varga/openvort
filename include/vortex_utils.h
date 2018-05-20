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
		int Npoints);
void insert_random_loops(struct tangle_state *tangle, int N);
void clip_at_wall(struct tangle_state *tangle);

void save_tangle(const char *filename, struct tangle_state *tangle);
void load_tangle(const char *filename, struct tangle_state *tangle);

int check_integrity(const struct tangle_state *tangle);
int is_empty(const struct tangle_state *tangle, int k);

#endif
