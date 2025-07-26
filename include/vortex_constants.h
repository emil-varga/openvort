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

extern int global_max_steps;

extern int global_LIA_only;

extern double global_alpha;
extern double global_alpha_p;

extern double global_dt;
extern double global_dl_min;
extern double global_dl_max;
extern double global_reconnection_angle_cutoff;
extern double rec_dist;

extern int global_small_loop_cutoff;
extern int global_frame_shot; //how often to save a snapshot of the tangle
extern int global_num_threads;

extern int global_use_mutual_friction;

extern int global_pin_mode;

extern int global_eliminate_origin_loops;
extern double global_eliminate_loops_origin_cutoff;

extern int global_eliminate_zaxis_loops;
extern double global_eliminate_loops_zaxis_cutoff;

extern int global_eliminate_outer_loops;
extern double global_eliminate_outer_loops_cutoff;

//injecting loops at upper z-plane
extern int global_loop_injection; //bool, inject or not
extern double global_loop_injection_frequency; //injections per second

//line pair injection
extern int global_line_injection; //bool, inject or not
extern int global_line_injection_n; //how many pairs to inject
extern double global_line_injection_frequency; //injections per second
extern int global_line_injection_polarized; //whether the injection is polarized

//Barnes-Hut ('tree') approximation
extern int global_use_BH;
extern double global_BH_resolution;
extern int global_BH_quadtree;
extern double global_BH_grain;

//hyperfriction to damp out strongly curved parts that are not resolved
//by the discretisation
extern int global_hyperfriction;
extern double global_max_curvature_scale;
extern double global_hyperalpha;

#endif //VORTEX_CONSTANTS_H
