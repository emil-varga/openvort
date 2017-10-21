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
  int curr_point, last_point, first_point;
  struct vec3d zdir = perpendicular(dir);
  
  vec3_cross(&u, &zdir, dir);
  vec3_cross(&v, &u, dir);

  vec3_mul(&u, &u, 1/vec3_d(&u));
  vec3_mul(&v, &v, 1/vec3_d(&v));

  printf("u:(%f, %f, %f), v:(%f, %f, %f)\n",
	 u.p[0], u.p[1], u.p[2],
	 v.p[0], v.p[1], v.p[2]);

  first_point = -1;

  for(int k=0; k<Npoints; ++k)
    {
      double phi = 2*M_PI/Npoints * k;
      double c = cos(phi);
      double s = sin(phi);
      vec3_mul(&p, &u, r*c);
      vec3_mul(&ptmp, &v, r*s);

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

void write_vector(FILE *stream, struct vec3d *v)
{
  fprintf(stream, "%.15g\t%.15g\t%.15g",
	  v->p[0], v->p[1], v->p[2]);
}

void save_point(FILE *stream, int vort_idx,
		const struct tangle_state *tangle, int i)
{
  fprintf(stream, "%d\t", vort_idx);
  write_vector(stream, tangle->vnodes + i);
  fprintf(stream, "\t");
  write_vector(stream, tangle->vels + i);
  fprintf(stream, "\t");
  write_vector(stream, tangle->tangents + i);
  fprintf(stream, "\t");
  write_vector(stream, tangle->normals + i);
  fprintf(stream, "\n");
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
	      save_point(stream, vortex_idx, tangle, curr);
	      visited[curr] = 1;
	      curr = tangle->connections[curr].forward;
	    }
	  save_point(stream, vortex_idx, tangle, curr);
	  visited[curr] = 1;
	  vortex_idx++;
	}
    }

  free(visited);
  fclose(stream);
}
