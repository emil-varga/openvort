#ifndef VORTEX_CONSTANTS_H
#define VORTEX_CONSTANTS_H

//these are defined in vortex_constants.c

extern double VORTEX_WIDTH;
extern double KAPPA;

extern double alpha;
extern double alpha_p;

/*
 * box sizes in the three directions
 *
 * these only apply for periodic/wall boundary conditions
 * the array is [xL, xH, yL, yH, zL, zH]
 */
extern double computation_box[6];

extern int use_mutual_friction;

#endif //VORTEX_CONSTANTS_H
