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

#include "configuration.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <math.h>

#include <libconfig.h>

#include "vec3_maths.h"
#include "vortex_utils.h"
#include "util.h"
#include "external_velocity.h"

#define PATH_LEN 256
//+1 for the \0
char output_dir[PATH_LEN+1];
char conf_file[PATH_LEN+1];

struct option options[] = {
    {"conf", required_argument, NULL, 'c'},
    {"output", required_argument, NULL, 'o'},
    {0, 0, 0, 0}
};

void set_walls_full(struct tangle_state *tangle, wall_type wall)
{
  for(int k=0; k<6; ++k)
        tangle->box.wall[k] = wall;
}

void set_walls_xy_periodic(struct tangle_state *tangle)
{
  set_walls_full(tangle, WALL_PERIODIC);
  tangle->box.wall[Z_L] = tangle->box.wall[Z_H] = WALL_OPEN;
}

int load_conf_vector(config_t *cfg, const char *name, struct vec3d *v)
{
  config_setting_t *vector;
  vector = config_lookup(cfg, name);
  if(vector && config_setting_type(vector) == CONFIG_TYPE_LIST)
    {
	*v = vec3(
	    config_setting_get_float_elem(vector, 0),
	    config_setting_get_float_elem(vector, 1),
	    config_setting_get_float_elem(vector, 2)
	  );
    }
  else
    goto failure;

  return 1;

failure:
  return 0;
}

int parse_options(int argc, char **argv)
{
  int c;
  int have_outdir = 0;
  int have_conf = 0;

  while(1)
    {
      c = getopt_long(argc, argv, "c:o:", options, NULL);

      if(c==-1) //all options parsed
	break;

      switch(c)
      {
	case 'c':
	  have_conf = 1;
	  strncpy(conf_file, optarg, PATH_LEN);
	  break;
	case 'o':
	  have_outdir = 1;
	  strncpy(output_dir, optarg, PATH_LEN);
	  break;
	default:
	  goto failure;
	  break;
      }
    }

  if(!have_outdir || !have_conf)
    goto failure;

  //trim the trailing / off of output_dir, if it's there
  char *x = output_dir;
  while(*(x+1)) x++;
  if(*x == '/') *x = 0;

  return 1;

failure:
  print_usage(argv[0]);
  return 0;
}

int setup_external_velocity(config_setting_t *v_conf_setting, struct v_conf_t *v_conf)
{
  double dval;
  const char *str;
  struct vec3d vval;
  config_setting_t *vec_cfg;

  config_setting_lookup_string(v_conf_setting, "type", &str);

  struct v_conf_t *test = &v_confs[0];
  for(;strlen(test->name) > 0; test++)
    {
      if(strcmp(test->name, str) == 0)
	{
	  printf("found v_conf: %s %s\n", test->name, str);
	  *v_conf = *test;
	  break;
	}
    }
  if(strlen(test->name) == 0)
    return -1; //we did not find the requested config

  //we found the requested config, now load its parameters, if any
  //TODO: error checking
  for(int k=0; k<(v_conf)->n_params; ++k)
    {
      switch((v_conf)->v_params[k].type)
      {
	case scalar_param:
	  config_setting_lookup_float(v_conf_setting, (v_conf)->v_params[k].name, &dval);
	  (v_conf)->v_params[k].value.scalar = dval;
	  break;
	case vector_param:
	  vec_cfg = config_setting_lookup(v_conf_setting, (v_conf)->v_params[k].name);
	  for(int j = 0; j<3; ++j)
	    vval.p[j] = config_setting_get_float_elem(vec_cfg, j);
	  (v_conf)->v_params[k].value.vector = vval;
	  break;
	default:
	  break;
      }
    }
  return 0;
}

int load_conf(const char *conf_file, struct tangle_state *tangle)
{
  config_t cfg;

  config_init(&cfg);
  if(!config_read_file(&cfg, conf_file))
    {
        fprintf(stderr, "%s:%d - %s\n", config_error_file(&cfg),
                config_error_line(&cfg), config_error_text(&cfg));
        config_destroy(&cfg);
        return 0;
    }

  int ival;
  double dval;
  const char *str;

  //setup general properties

  if(config_lookup_bool(&cfg, "use_mutual_friction", &ival))
    use_mutual_friction = ival;

  if(config_lookup_float(&cfg, "alpha", &dval))
    alpha = dval;

  if(config_lookup_float(&cfg, "alpha_p", &dval))
    alpha_p = dval;

  if(config_lookup_float(&cfg, "KAPPA", &dval))
    KAPPA = dval;

  if(config_lookup_float(&cfg, "dt", &dval))
    global_dt = dval;

  if(config_lookup_float(&cfg, "dl_min", &dval))
    global_dl_min = dval;
  if(config_lookup_float(&cfg, "dl_max", &dval))
    global_dl_max = dval;
  if(config_lookup_int(&cfg, "small_loop_cutoff", &ival))
    small_loop_cutoff = ival;
  if(config_lookup_float(&cfg, "reconnection_angle_cutoff", &dval))
    reconnection_angle_cutoff = M_PI/180.0 * dval;
  if(config_lookup_float(&cfg, "reconnection_distance", &dval))
     rec_dist = dval;
  if(config_lookup_int(&cfg, "frame_shots", &ival))
    frame_shot = ival;

  if(config_lookup_bool(&cfg, "eliminate_origin_loops", &ival))
      eliminate_origin_loops = ival;
  if(config_lookup_float(&cfg, "eliminate_loops_origin_cutoff", &dval))
      eliminate_loops_origin_cutoff = dval;

  if(config_lookup_bool(&cfg, "eliminate_zaxis_loops", &ival))
    eliminate_zaxis_loops = ival;
  if(config_lookup_float(&cfg, "eliminate_loops_origin_cutoff", &dval))
    eliminate_loops_zaxis_cutoff = dval;

  if(config_lookup_int(&cfg, "num_threads", &ival))
    global_num_threads = ival;

  config_setting_t *vel_conf;
  vel_conf = config_lookup(&cfg, "vn_conf");
  if(vel_conf && config_setting_type(vel_conf) == CONFIG_TYPE_GROUP)
    setup_external_velocity(vel_conf, &vn_conf);
  printf("set up %s\n", vn_conf.name);


  vel_conf = config_lookup(&cfg, "vs_conf");
  if(vel_conf && config_setting_type(vel_conf) == CONFIG_TYPE_GROUP)
    setup_external_velocity(vel_conf, &vs_conf);

  //setup the domain box size

  config_setting_t *domain, *point;
  domain = config_lookup(&cfg, "domain");
  if(domain && config_setting_type(domain) == CONFIG_TYPE_LIST)
    {
	  point = config_setting_get_elem(domain, 0);
	  tangle->box.bottom_left_back = vec3(
	      config_setting_get_float_elem(point, 0),
	      config_setting_get_float_elem(point, 1),
	      config_setting_get_float_elem(point, 2)
	      );

	  point = config_setting_get_elem(domain, 1);
	  tangle->box.top_right_front = vec3(
	      config_setting_get_float_elem(point, 0),
	      config_setting_get_float_elem(point, 1),
	      config_setting_get_float_elem(point, 2)
	      );
    }
  else
    {
      fprintf(stderr, "Error: incorrectly specified domain\n");
      goto failure;
    }

  //setup the boundaries

  if(!config_lookup_string(&cfg, "boundaries", &str))
    {
      fprintf(stderr, "Error: specify boundaries\n");
      goto failure;
    }
  if(strcmp(str, "periodic-6") == 0)
    {
      tangle->bimg = periodic_6;
      set_walls_full(tangle, WALL_PERIODIC);
    }
  else if(strcmp(str, "periodic-18") == 0)
    {
      tangle->bimg = periodic_18;
      set_walls_full(tangle, WALL_PERIODIC);
    }
  else if(strcmp(str, "periodic-26") == 0)
    {
      tangle->bimg = periodic_26;
      set_walls_full(tangle, WALL_PERIODIC);
    }
  else if(strcmp(str, "wall-1-6") == 0)
    {
      tangle->bimg = wall_1_6;
      set_walls_full(tangle, WALL_PERIODIC);
      tangle->box.wall[Z_L] = WALL_MIRROR;
      tangle->box.wall[Z_H] = WALL_OPEN;
    }
  else if(strcmp(str, "wall-1-18") == 0)
    {
      tangle->bimg = wall_1_18;
      set_walls_full(tangle, WALL_PERIODIC);
      tangle->box.wall[Z_L] = WALL_MIRROR;
      tangle->box.wall[Z_H] = WALL_OPEN;
    }
  else if(strcmp(str, "wall-1-26") == 0)
    {
      tangle->bimg = wall_1_26;
      set_walls_full(tangle, WALL_PERIODIC);
      tangle->box.wall[Z_L] = WALL_MIRROR;
      tangle->box.wall[Z_H] = WALL_OPEN;
    }
  else if(strcmp(str, "wall-1-open") == 0)
    {
      tangle->bimg = wall_1_open;
      set_walls_full(tangle, WALL_OPEN);
      tangle->box.wall[Z_L] = WALL_MIRROR;
    }
  else if(strcmp(str, "open") == 0)
    {
      tangle->bimg = open_boundaries;
      set_walls_full(tangle, WALL_OPEN);
    }
  else
    {
      fprintf(stderr, "Error: unknown boundary condition\n");
      goto failure;
    }

  //setup the initial condition

  if(!setup_init(conf_file, tangle))
    goto failure;

  config_destroy(&cfg);
  return 1;
failure:
  config_destroy(&cfg);
  return 0;
}

int setup_init(const char *conf_file, struct tangle_state *tangle)
{
  config_t cfg;

  config_init(&cfg);
  if(!config_read_file(&cfg, conf_file))
    {
      fprintf(stderr, "Can't read config file %s", conf_file);
      return 0;
    }

  const char *str;
  int ival;

  if(config_lookup_string(&cfg, "init_mode", &str))
    {
      if(strcmp(str, "random loops") == 0)
	{
	  if(config_lookup_int(&cfg, "init_loops_n", &ival))
	    insert_random_loops(tangle, ival);
	  else
	    {
	      fprintf(stderr, "Error: set the number of loops in config_file\n");
	      goto failure;
	    }
	  if(tangle->box.wall[Z_L] == WALL_MIRROR)
	    clip_at_wall(tangle);
	}
      else if(strcmp(str, "one loop") == 0)
	{
	  struct vec3d c;
	  struct vec3d d;
	  double r;
	  int N;
	  load_conf_vector(&cfg, "loop_center", &c);
	  load_conf_vector(&cfg, "loop_dir", &d);
	  if(!config_lookup_float(&cfg, "loop_r", &r))
	    goto failure;
	  if(!config_lookup_int(&cfg, "loop_N", &N))
	    goto failure;
	  add_circle(tangle, &c, &d, r, N);

	  if(tangle->box.wall[Z_L] == WALL_MIRROR)
	    clip_at_wall(tangle);
	}//TODO: add more init modes
      else
	{
	  fprintf(stderr, "Error: unknown initialization mode: %s\n", str);
	  goto failure;
	}
    }
  else
    {
      fprintf(stderr, "Error: set the initialization mode\n");
      goto failure;
    }

  config_destroy(&cfg);
  return 1;

failure:
  config_destroy(&cfg);
  return 0;
}

void print_ev_config(struct v_conf_t *v_conf)
{
  printf("type = %s\n", (v_conf)->name);
  printf("With parameters:\n");
  for(int k=0; k < (v_conf)->n_params; ++k)
    {
      printf("%s = ", (v_conf)->v_params[k].name);
      switch((v_conf)->v_params[k].type)
	{
	case scalar_param:
  	  printf("%g\n", (v_conf)->v_params[k].value.scalar);
  	  break;
  	case vector_param:
  	  printf("(%g, %g, %g)\n", (v_conf)->v_params[k].value.vector.p[0],
  		 (v_conf)->v_params[k].value.vector.p[1],
  		 (v_conf)->v_params[k].value.vector.p[2]);
  	  break;
        }
      }
}

void print_config()
{
  printf("Normal fluid setup:\n");
  print_ev_config(&vn_conf);

  printf("Superfluid setup:\n");
  print_ev_config(&vs_conf);

  printf(
      "time step                  = %g s\n"
      "minimum distance           = %g cm\n"
      "maximum distance           = %g cm\n"
      "minimum reconnection angle = %g rad\n"
      "alpha                      = %g\n"
      "alpha_p                    = %g\n"
      "small loop cutoff          = %d points\n"
      "steps per frame            = %d\n"
      "number of threads          = %d\n"
      "loop removal near origin   = %s\n"
      "loop removal cutoff        = %g cm\n",
      global_dt,
      global_dl_min,
      global_dl_max,
      reconnection_angle_cutoff,
      alpha,
      alpha_p,
      small_loop_cutoff,
      frame_shot,
      global_num_threads,
      eliminate_origin_loops ? "On" : "Off",
      eliminate_loops_origin_cutoff);
}
