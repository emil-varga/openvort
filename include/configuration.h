/*
 * configuration.h
 *
 *  Created on: Mar 8, 2018
 *      Author: emil
 *
 *
 *  The configuration support is handled through libconfig.
 *
 *  Functions declared here load the configuration file, modify the
 *  global constants declared in vortex_constants.h and set up the
 *  initial condition
 */

#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include "vortex_constants.h"
#include "tangle.h"

int load_conf(const char *conf_file, struct tangle_state *tangle);
int setup_init(const char *conf_file, struct tangle_state *tangle);

#endif /* CONFIGURATION_H */
