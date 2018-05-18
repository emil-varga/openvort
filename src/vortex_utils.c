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
    {
      return;
    }
  if(!Npoints)
    return;

  struct vec3d u, v, p, ptmp;
  int curr_point, last_point, first_point;
  struct vec3d zdir = perpendicular(dir);
  
  vec3_cross(&u, &zdir, dir);
  vec3_cross(&v, dir, &u);

  vec3_mul(&u, &u, 1/vec3_d(&u));
  vec3_mul(&v, &v, 1/vec3_d(&v));

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

void insert_random_loops(struct tangle_state *tangle, int N)
{
  //TODO: configurable range of loop radii?, number of points?
  double D = 1;
  double Ls[3];
  for(int k = 0; k<3; ++k)
    {
	Ls[k] = tangle->box.top_right_front.p[k] - tangle->box.bottom_left_back.p[k];
	if(Ls[k] < D)
	  D = Ls[k];
    }

  const double rmin = 0.05*D;
  const double rmax = D;

  for(int k = 0; k<N; ++k)
    {
      struct vec3d dir = vec3(drand48()-0.5, drand48()-0.5, drand48()-0.5);
      vec3_normalize(&dir);

      struct vec3d c;
      for(int j=0; j<3; ++j)
	c.p[j] = tangle->box.bottom_left_back.p[j] + drand48()*Ls[j];

      double r =  rmin + drand48()*(rmax - rmin);
      add_circle(tangle, &c, &dir, r, 64);
    }
}


void clip_at_wall(struct tangle_state *tangle)
{
  /*
   * Clips the tangle at the lower z-wall
   */
  double limit = tangle->box.bottom_left_back.p[2];
  for(int k=0; k<tangle->N; ++k)
    tangle->recalculate[k] = 0;

  for(int kk=0; kk<tangle->N; ++kk)
    {
      if(tangle->status[kk].status == EMPTY ||
	  tangle->recalculate[kk] > 0)
	continue;

      int here = kk;
      int next = tangle->connections[kk].forward;
      int tmp;

#define lchk(z) (tangle->vnodes[z].p[2] > limit)
      int seek_under = lchk(here);

      while(next != kk)
	{
	  if(seek_under) //we are above and are looking for the crossing point
	    {
	      if(!lchk(next)) //the next is below
		{
		  seek_under = !seek_under;
		  tangle->vnodes[next].p[2] = limit;
		  tangle->recalculate[here]++;
		  tangle->recalculate[next]++;
		  tmp = tangle->connections[next].forward;

		  tangle->connections[next].forward = -1;
		  tangle->status[next].status = PINNED;
		  tangle->status[next].pin_wall = Z_L;

		  here = tmp;
		  next = tangle->connections[here].forward;
		  continue;
		}
	      tangle->recalculate[here]++;
	      tangle->recalculate[next]++;
	      here = next;
	      next = tangle->connections[here].forward;
	    }
	  else //we are bellow and looking for the crossing point
	    {
	      if(lchk(next)) //the next is above
		{
		  seek_under = !seek_under;
		  tangle->vnodes[here].p[2] = limit;
		  tangle->recalculate[here]++;
		  tangle->recalculate[next]++;
		  tangle->connections[here].reverse = -1;
		  tangle->status[here].status = PINNED;
		  tangle->status[here].pin_wall = Z_L;
		  here = next;
		  next = tangle->connections[here].forward;
		  continue;
		}
	      tangle->recalculate[here]++;
	      tangle->recalculate[next]++;
	      tangle->status[here].status = EMPTY;
	      tangle->status[next].status = EMPTY;
	      here = next;
	      next = tangle->connections[here].forward;
	    }
	}
    }
#undef lchk
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

int load_point(FILE *stream, struct tangle_state *tangle, int idx)
{
  int vidx;
  struct vec3d *pos     = &tangle->vnodes[idx];
  struct vec3d *vel     = &tangle->vels[idx];
  struct vec3d *tangent = &tangle->tangents[idx];
  struct vec3d *normal  = &tangle->normals[idx];
  int ret = fscanf(stream, "%d\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg",
	   &vidx,
	   pos->p,     pos->p+1,     pos->p+2,
	   vel->p,     vel->p+1,     vel->p+2,
	   tangent->p, tangent->p+1, tangent->p+2,
	   normal->p,  normal->p+1,  normal->p+2);
  if(ret == EOF)
    return -1;

  tangle->status[idx].status = FREE;

  return vidx;
}

void save_tangle(const char *filename, struct tangle_state *tangle)
{
  int vortex_idx = 0;
  int *visited = (int *)calloc(tangle->N, sizeof(int));
  FILE *stream = fopen(filename, "w");

  if(!stream) {
      fprintf(stderr, "Can't open file %s\n", filename);
      return;
  }

  for(int k=0; k < tangle->N; ++k)
    {
      if(!visited[k])
	{
	  if(tangle->status[k].status == EMPTY)
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
	      if(curr < 0)
		{
		  curr = tangle->connections[first].reverse;
		  while(tangle->connections[curr].reverse > 0)
		    {
		      save_point(stream, vortex_idx, tangle, curr);
		      visited[curr] = 1;
		      curr = tangle->connections[curr].reverse;
		    }
		  break;
		}
	    }
	  if(curr > 0)
	    {
	      save_point(stream, vortex_idx, tangle, curr);
	      visited[curr] = 1;
	    }
	  vortex_idx++;
	}
    }

  free(visited);
  fclose(stream);
}

/*
 * Loads a tangle from a file.
 * TODO: does not support walls yet, only closed loops
 */
int load_tangle(const char *filename, struct tangle_state *tangle)
{
  FILE *file = fopen(filename, "r");
  int vidx = 0; //vortex index
  int current_vidx; //index of vortex we are loading now
  int pidx = 0; //index of position to store a loaded point
  int pidx_start = 0; //where does the vortex begin
  while((current_vidx = load_point(file, tangle, pidx)) >= 0)
    {
      tangle->connections[pidx].reverse = pidx - 1;
      tangle->connections[pidx].forward = pidx + 1;

      if(current_vidx != vidx)
	{
	  //we are on a new vortex
	  //stitch together the vortex we finished
	  tangle->connections[pidx_start].reverse = pidx-1;
	  tangle->connections[pidx-1].forward = pidx_start;

	  //move on to the next
	  vidx++;
	  pidx_start = pidx;
	  if(current_vidx != vidx)
	    return -1; //something's wrong

	}

      pidx++;
      if(pidx > tangle->N - 1)
      	expand_tangle(tangle, 2*tangle->N);
    }

  return 0;
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
      if(next >= 0 && tangle->status[k].status == FREE)
	{
	  if(k != tangle->connections[next].reverse)
	    error("Forward connection broken %d %d %d\n", k, next,
		  tangle->connections[next].reverse);
	}
      if(prev >= 0 && tangle->status[k].status == FREE)
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
