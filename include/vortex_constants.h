#ifndef VORTEX_CONSTANTS_H
#define VORTEX_CONSTANTS_H

//these are defined with default values in vortex_constants.c

extern double VORTEX_WIDTH;
extern double KAPPA;

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

#endif //VORTEX_CONSTANTS_H
