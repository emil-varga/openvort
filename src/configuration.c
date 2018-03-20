/*
 * configuration.c
 *
 *  Created on: Mar 8, 2018
 *      Author: emil
 */

#include "configuration.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>

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

int setup_external_velocity(config_setting_t *v_conf_setting, struct v_conf_t **v_conf)
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
	  *v_conf = test;
	  break;
	}
    }
  if(strlen(test->name) == 0)
    return -1; //we did not find the requested config

  //we found the requested config, now load its parameters, if any
  //TODO: erro checking
  for(int k=0; k<(*v_conf)->n_params; ++k)
    {
      switch((*v_conf)->v_params[k].type)
      {
	case scalar_param:
	  config_setting_lookup_float(v_conf_setting, (*v_conf)->v_params[k].name, &dval);
	  (*v_conf)->v_params[k].value.scalar = dval;
	  break;
	case vector_param:
	  vec_cfg = config_setting_lookup(v_conf_setting, (*v_conf)->v_params[k].name);
	  for(int j = 0; j<3; ++j)
	    vval.p[j] = config_setting_get_float_elem(vec_cfg, j);
	  (*v_conf)->v_params[k].value.vector = vval;
	  break;
	default:
	  break;
      }
    }

  printf("Using external velocity %s\n", (*v_conf)->name);
  printf("With parameters:\n");
  for(int k=0; k < (*v_conf)->n_params; ++k)
    {
      printf("%s = ", (*v_conf)->v_params[k].name);
      switch((*v_conf)->v_params[k].type)
      {
	case scalar_param:
	  printf("%g\n", (*v_conf)->v_params[k].value.scalar);
	  break;
	case vector_param:
	  printf("(%g, %g, %g)\n", (*v_conf)->v_params[k].value.vector.p[0],
		 (*v_conf)->v_params[k].value.vector.p[1],
		 (*v_conf)->v_params[k].value.vector.p[2]);
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


  config_setting_t *vel_conf;
  vel_conf = config_lookup(&cfg, "vn_conf");
  if(vel_conf && config_setting_type(vel_conf) == CONFIG_TYPE_GROUP)
    setup_external_velocity(vel_conf, &vn_conf);

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
    }
  else if(strcmp(str, "wall-1-18") == 0)
    {
      tangle->bimg = wall_1_18;
      set_walls_full(tangle, WALL_PERIODIC);
    }
  else if(strcmp(str, "wall-1-26") == 0)
    {
      tangle->bimg = wall_1_26;
      set_walls_full(tangle, WALL_PERIODIC);
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
	} //TODO: add more init modes
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

int setup_external_velocity(config_setting_t *v_conf_setting, struct v_conf_t **v_conf);
