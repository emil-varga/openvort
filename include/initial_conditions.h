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

#endif /* INCLUDE_INITIAL_CONDITIONS_H_ */
