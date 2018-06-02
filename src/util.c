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

#include "util.h"
#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <vec3_maths.h>

void error(const char *msg, ...)
{
  va_list args;
  va_start(args, msg);
  vfprintf(stderr, msg, args);
  abort();
}

void print_vec(const struct vec3d *v)
{
  printf("(%.3g, %.3g, %.3g)", v->p[0], v->p[1], v->p[2]);
}

void print_usage(const char *prog_name)
{
  printf("Usage:\n");
  fputs(prog_name, stdout);
  printf(" -c | --conf <path to config file> -o | --output <output directory>\n");
}

int setup_outdir(const char *dirname)
{
  DIR *test = opendir(dirname);

  if(test){
      //the directory exist and we can open it
      closedir(test);
      return 1;
  } else {
      //either it doesn't exist or we can't open it
      //try creating it
      if(mkdir(dirname, 0777) == -1)
	return 0;
  }

  return 1;
}
