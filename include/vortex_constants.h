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

#ifndef VORTEX_CONSTANTS_H
#define VORTEX_CONSTANTS_H

//these are defined with default values in vortex_constants.c

extern double VORTEX_WIDTH;
extern double KAPPA;

extern int max_steps;

extern int LIA_only;

extern double alpha;
extern double alpha_p;

extern double global_dt;
extern double global_dl_min;
extern double global_dl_max;
extern double reconnection_angle_cutoff;
extern double rec_dist;

extern int small_loop_cutoff;
extern int frame_shot; //how often to save a snapshot of the tangle
extern int global_num_threads;

extern int use_mutual_friction;

extern int pin_mode;

extern int eliminate_origin_loops;
extern double eliminate_loops_origin_cutoff;

extern int eliminate_zaxis_loops;
extern double eliminate_loops_zaxis_cutoff;

extern int eliminate_outer_loops;
extern double eliminate_outer_loops_cutoff;

//injecting loops at upper z-plane
extern int loop_injection; //bool, inject or not
extern double loop_injection_frequency; //injections per second

//line pair injection
extern int line_injection; //bool, inject or not
extern int line_injection_n; //how many pairs to inject
extern double line_injection_frequency; //injections per second
extern int line_injection_polarized; //whether the injection is polarized

//Barnes-Hut ('tree') approximation
extern int use_BH;
extern double BH_resolution;
extern int BH_quadtree;
extern double BH_grain;

//hyperfriction to damp out strongly curved parts that are not resolved
//by the discretisation
extern int hyperfriction;
extern double max_curvature_scale;
extern double hyperalpha;

#endif //VORTEX_CONSTANTS_H
