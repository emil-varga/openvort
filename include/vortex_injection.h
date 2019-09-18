/*
 * vortex_injection.h
 *
 *  Created on: Sep 10, 2019
 *      Author: emil
 */

#ifndef INCLUDE_VORTEX_INJECTION_H_
#define INCLUDE_VORTEX_INJECTION_H_

#include "tangle.h"

void inject_vortices(struct tangle_state *tangle, double t);
void inject_loop(struct tangle_state *tangle, double t, double frequency);
void inject_line_pairs(struct tangle_state *tangle, double t, double frequency, int line_injection_n,
		       int polarized);


#endif /* INCLUDE_VORTEX_INJECTION_H_ */
