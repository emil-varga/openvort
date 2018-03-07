#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "vortex_utils.h"
#include "util.h"

struct vec3d perpendicular(const struct vec3d *dir)
{
  struct vec3d a = vec3(1,0,0);
  struct vec3d b = vec3(0,1,0);
  struct vec3d res;

  vec3_mul(&res, dir, -vec3_dot(&a, dir));
  vec3_add(&res, &res, &a);

  if(vec3_d(&res) > 1e-8)
    return res;

  vec3_mul(&res, dir, -vec3_dot(&b, dir));
  vec3_add(&res, &res, &b);
  return res;
}  

void add_circle(struct tangle_state *tangle,
		struct vec3d *center, struct vec3d *dir, double r,
		int Npoints)
{
  if(num_free_points(tangle) < Npoints)
    return;
  if(!Npoints)
    return;

  struct vec3d u, v, p, ptmp;
  int curr_point, last_point, first_point;
  struct vec3d zdir = perpendicular(dir);
  
  vec3_cross(&u, &zdir, dir);
  vec3_cross(&v, dir, &u);

  vec3_mul(&u, &u, 1/vec3_d(&u));
  vec3_mul(&v, &v, 1/vec3_d(&v));

  printf("u:(%f, %f, %f), v:(%f, %f, %f)\n",
	 u.p[0], u.p[1], u.p[2],
	 v.p[0], v.p[1], v.p[2]);

  first_point = -1;
  curr_point = -1;

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

      tangle->status[curr_point].status = FREE;
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

int check_loop(const struct tangle_state *tangle, int *visited, int k);
int check_integrity(const struct tangle_state *tangle)
{
  int *visited = calloc(tangle->N, sizeof(int));
  int errors = 0;

  for(int k=0; k < tangle->N; ++k)
    {
      int next = tangle->connections[k].forward;
      int prev = tangle->connections[k].reverse;
      if(next >= 0)
	{
	  if(k != tangle->connections[next].reverse)
	    error("Forward connection broken %d %d %d\n", k, next,
		  tangle->connections[next].reverse);
	}
      if(prev >= 0)
	{
	  if( k!= tangle->connections[prev].forward)
	    error("Reverse connection broken %d %d %d\n", k, prev,
		  tangle->connections[prev].forward);
	}
      if(!visited[k])
	{
	  visited[k] += 1;	  	  
	  errors += check_loop(tangle, visited, k);
	}
    }

  free(visited);
  return errors;
}

int check_loop(const struct tangle_state *tangle, int *visited, int k)
{
  int errors = 0;
  int rval;

  if((rval = is_empty(tangle, k)) > 0)
    return 0;
  if(rval < 0)
    return rval;
  
  int j = tangle->connections[k].forward;
  int total = 0;
  while(j != k && total < tangle->N)
    {
      if(visited[j])
	{
#ifdef _DEBUG_
	  printf("Ran into point %d again.\n", j);
#endif
	  errors--;
	  break;
	}
      visited[j]++;
      if((rval = is_empty(tangle, j)) > 0)
	{
#ifdef _DEBUG_
	  printf("Connected to an empty point, %d %d.\n",
		 k, j);
#endif
	  errors--;
	  return errors;
	}
      if(tangle->connections[j].forward < 0)
	{
#ifdef _DEBUG_
	  printf("Linked point with -1 forward, %d.\n",
		 j);
#endif
	  errors--;
	  return errors;
	}

      j = tangle->connections[j].forward;

      total++;
    }


  if(j != k && total == tangle->N)
    {
#ifdef _DEBUG_
      printf("Ran through the entire tangle from %d.\n", k);
#endif
      errors--;
    }

  return errors;
}

int is_empty(const struct tangle_state *tangle, int k)
{
  struct neighbour_t nb = tangle->connections[k];

  if(nb.forward == -1 &&
     nb.reverse == -1)
    return 1; //point is empty

  if(nb.forward >= 0 &&
     nb.reverse >= 0)
    return 0; //both links are valid, point is not empty

  //one of the links is >= 0, the other is ==-1
  //point is corrupted
  #ifdef _DEBUG_
  printf("Corrupted links in point %d (%d, %d)\n",
	 k,
	 tangle->connections[k].forward,
	 tangle->connections[k].reverse);
  #endif
  return -1;
}
