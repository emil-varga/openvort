#include <math.h>
#include "vortex_utils.h"
#include <stdio.h>

struct vec3d perpendicular(const struct vec3d *dir)
{
  struct vec3d a = vec3(1,0,0);
  struct vec3d b = vec3(0,1,0);
  struct vec3d res;

  vec3_mul(&res, dir, -vec3_dot(&a, dir));
  vec3_add(&res, &res, &a);

  if(vec3_d(&res) > 0)
    return res;

  vec3_mul(&res, dir, -vec3_dot(&b, dir));
  vec3_add(&res, &res, &b);
  return res;
}  

void add_circle(struct tangle_state *tangle,
		struct vec3d *center, struct vec3d *dir, double r,
		size_t Npoints)
{
  if(num_free_points(tangle) < Npoints)
    return;

  struct vec3d u, v, p, ptmp;
  int curr_point, last_point, first_point, k;
  struct vec3d zdir = perpendicular(dir);
  
  vec3_cross(&u, &zdir, dir);
  vec3_cross(&v, &u, dir);

  vec3_mul(&u, &u, 1/vec3_d(&u));
  vec3_mul(&v, &v, 1/vec3_d(&v));

  printf("u:(%f, %f, %f), v:(%f, %f, %f)\n",
	 u.p[0], u.p[1], u.p[2],
	 v.p[0], v.p[1], v.p[2]);

  first_point = -1;

  for(k=0; k<Npoints; ++k)
    {
      double phi = 2*M_PI/Npoints * k;
      vec3_mul(&p, &u, r*cos(phi));
      vec3_mul(&ptmp, &v, r*sin(phi));

      vec3_add(&p, &p, &ptmp);
      vec3_add(&p, &p, center);
      curr_point = get_tangle_next_free(tangle);
      if(first_point < 0)
	first_point = curr_point;

      tangle->vnodes[curr_point] = p;
      if(curr_point != first_point)
	{
	  tangle->connections[curr_point].reverse = last_point;
	  tangle->connections[last_point].forward = curr_point;
	}

      last_point = curr_point;
    }
  tangle->connections[curr_point].forward  = first_point;
  tangle->connections[first_point].reverse = curr_point;
}

void save_tangle(const char *filename, struct tangle_state *tangle)
{
  int vortex_idx = 0;
  int *visited = (int *)calloc(tangle->N, sizeof(int));
  FILE *stream = fopen(filename, "w");

  for(int k=0; k < tangle->N; ++k)
    {
      if(!visited[k])
	{
	  if(tangle->connections[k].forward < 0)
	    {
	      visited[k] = 1;
	      continue;
	    }
	  
	  int first = k;
	  int curr  = k;
	  while(tangle->connections[curr].forward != first)
	    {
	      fprintf(stream, "%d\t%g\t%g\t%g\n", vortex_idx,
		      tangle->vnodes[curr].p[0],
		      tangle->vnodes[curr].p[1],
		      tangle->vnodes[curr].p[2]);
	      visited[curr] = 1;
	      curr = tangle->connections[curr].forward;
	    }
	  fprintf(stream, "%d\t%g\t%g\t%g\n", vortex_idx,
		  tangle->vnodes[curr].p[0],
		  tangle->vnodes[curr].p[1],
		  tangle->vnodes[curr].p[2]);
	  visited[curr] = 1;
	  vortex_idx++;
	}
    }

  free(visited);
  fclose(stream);
}
