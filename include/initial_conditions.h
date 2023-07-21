/*
 * initial_conditions.h
 *
 *  Created on: Sep 10, 2019
 *      Author: emil
 */

#ifndef INCLUDE_INITIAL_CONDITIONS_H_
#define INCLUDE_INITIAL_CONDITIONS_H_

#include "tangle.h"

void random_straight_lines(struct tangle_state *tangle, int npairs, int points_per_line);
void insert_random_loops(struct tangle_state *tangle, int N);
void make_big_ring(struct tangle_state *tangle, double ring_r, int ring_N);
void quad_straight_lines(struct  tangle_state *tangle, int points_per_line);
void dipole_straight_lines(struct  tangle_state *tangle, int points_per_line, int direction);
void single_straight_line(struct tangle_state *tangle, int points_per_line, int direction);
void random_polarized_lines(struct tangle_state *tangle, int nvort, int direction, int points_per_line);

#endif /* INCLUDE_INITIAL_CONDITIONS_H_ */
