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

#ifndef VORTEX_DYNAMICS
#define VORTEX_DYNAMICS

#include "tangle.h"

/*
  All functions dealing with movement of the vortices should go here.
  
  Particularly, stepping and reconnecting.
 */


// Stepping
struct rk4_state {
  //the tangle states here should be reduced
  //only calculate what we really need for the velocity
  struct tangle_state *s1, *s2, *s3, *s4;

  double dt;
};

struct ab4_state {
  struct tangle_state *s1, *s2, *s3, *s4;

  double dt;
};

void euler_step(struct tangle_state *tangle, double dt);

void rk4_step(struct tangle_state *tangle, double dt);
void ab4_init_state(struct ab4_state **state,
		    const struct tangle_state *tangle, double dt);
void ab4_step(struct tangle_state *tangle, struct ab4_state *state);


//reconnections

int reconnect(struct tangle_state *tangle, double rec_dist, double rec_angle);

#endif //VORTEX_DYNAMICS
