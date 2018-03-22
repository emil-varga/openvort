/*
 * vortex_constants.c
 *
 *  Created on: Dec 25, 2017
 *      Author: emil
 */

#include "vortex_constants.h"
//declared as extern in vortex_constants.h

double VORTEX_WIDTH = 1e-8; //cm
//quantum of circulation
double KAPPA = 9.97e-4; //cm^2/s

//mutual friction parameters
double alpha = 0.1;
double alpha_p = 0.01;

//densities, g/cm^3
double rho_n = 1.0; //normal fluid
double rho_s = 0.45; //superfluid

//discretisation and stepping
double global_dt = 1e-3;
double global_dl_min = 1e-3;
double global_dl_max = 5e-3;
double reconnection_angle_cutoff = 0.087; //in radians, about 5 deg
double rec_dist = 1e-3;

int small_loop_cutoff = 5;
int frame_shot = 100;
int global_num_threads = 4;

//switch mutual friction on/off
int use_mutual_friction = 1;
