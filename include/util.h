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

#ifndef UTIL_H
#define UTIL_H

#include <assert.h>
#include <math.h>

/*
 * Assertion that also prints a useful message.
 */
#ifndef NDEBUG
#define assert_msg(expr,msg)\
  if(!(expr)) {\
      printf(msg);\
      assert((expr));\
  }
#else
#define assert_msg(expr,msg) {}
#endif //ifndef NDEBUG

void error(const char *msg, ...);
void print_usage(const char *prog_name);
int setup_outdir(const char *dirname);

#define DEG2RAD(X) ((X)*M_PI/180.0)

#endif
