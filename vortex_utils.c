#include <math.h>
#include "vortex_utils.h"

void add_circle(struct tangle_state *tangle,
		struct vec3d *center, struct vec3d *dir, double r,
		size_t Npoints)
{
  if(tangle->N - tangle->next_free < Npoints)
    return;

  struct vec3d u, v, p, ptmp;
  int curr_point, last_point, first_point;
  vec3_cross(&u, center, dir);
  vec3_cross(&v, &u, dir);

  vec3_mul(&u, &u, 1/vec3_d(&u));
  vec3_mul(&v, &v, 1/vec3_d(&v));

  first_point = last_point = curr_point  = tangle->next_free;

  for(double phi = 0; phi < 2*M_PI; phi += 2*M_PI/Npoints)
    {
      vec3_mul(&p, &u, r*cos(phi));
      vec3_mul(&ptmp, &v, r*sin(phi));

      vec3_add(&p, &p, &ptmp);
      curr_point = tangle->next_free;

      tangle->vnodes[curr_point] = p;
      if(curr_point != first_point)
	{
	  tangle->connections[curr_point].reverse = last_point;
	  tangle->connections[last_point].forward = curr_point;
	}

      last_point = curr_point;
      update_tangle_next_free(tangle);
    }
  tangle->connections[curr_point].forward  = first_point;
  tangle->connections[first_point].reverse = curr_point;
}

void save_tangle(FILE *stream, struct tangle_state *tangle)
{
  int vortex_idx = 0;

  int *visited = (int *)calloc(tangle->N, sizeof(int))
}
