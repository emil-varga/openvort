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
#include "normal_fluid.h"

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

  if(!config_lookup_string(&cfg, "normal-vel", &str))
    get_vn = get_vn_noflow;
  else
    {
      if(strcmp(str, "no flow") == 0)
	get_vn = get_vn_noflow;
      else if(strcmp(str, "spherical") == 0)
	{
	  get_vn = get_vn_spherical;
	  config_lookup_float(&cfg, "spherical_vn_strength", &dval);
	  spherical_vn_strength = dval;
	  config_lookup_float(&cfg, "spherical_vn_cutoff", &dval);
	  spherical_vn_cutoff = dval;
	}
    }

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

void insert_random_loops(struct tangle_state *tangle, int N);


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

void insert_random_loops(struct tangle_state *tangle, int N)
{
  //TODO: configurable range of loop radii?, number of points?
  double Ls[3];
  for(int k = 0; k<3; ++k)
	Ls[k] = tangle->box.top_right_front.p[k] - tangle->box.bottom_left_back.p[k];

  const double D = (Ls[0] + Ls[1] + Ls[2]) / 3.0;
  const double rmin = 0.05*D;
  const double rmax = D;

  for(int k = 0; k<N; ++k)
    {
      struct vec3d dir = vec3(drand48(), drand48(), drand48());
      vec3_normalize(&dir);

      struct vec3d c;
      for(int j=0; j<3; ++j)
	c.p[j] = tangle->box.bottom_left_back.p[j] + drand48()*Ls[j];

      double r =  rmin + drand48()*(rmax - rmin);
      add_circle(tangle, &c, &dir, r, 64);
    }
}
