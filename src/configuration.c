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

#include <libconfig.h>

#include "vec3_maths.h"
#include "vortex_utils.h"

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
    tangle->bimg = periodic_6;
  else if(strcmp(str, "periodic-18") == 0)
    tangle->bimg = periodic_18;
  else if(strcmp(str, "periodic-26") == 0)
    tangle->bimg = periodic_26;
  else if(strcmp(str, "wall-1-6") == 0)
    tangle->bimg = wall_1_6;
  else if(strcmp(str, "wall-1-18") == 0)
    tangle->bimg = wall_1_18;
  else if(strcmp(str, "wall-1-26") == 0)
    tangle->bimg = wall_1_26;
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
