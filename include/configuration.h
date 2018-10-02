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
 ******************************************************************************

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

extern char output_dir[];
extern char conf_file[];

int parse_options(int argc, char **argv);
int load_conf(const char *conf_file, struct tangle_state *tangle);
int setup_init(const char *conf_file, struct tangle_state *tangle);
void print_config(const struct tangle_state *tangle);

#endif /* CONFIGURATION_H */
