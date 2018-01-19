/*
 * vortex_constants.c
 *
 *  Created on: Dec 25, 2017
 *      Author: emil
 */

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

//switch mutual friction on/off
int use_mutual_friction = 1;
