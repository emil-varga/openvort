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
#include <math.h>

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
    {
      vec3_normalize(&res);
      return res;
    }

  vec3_mul(&res, dir, -vec3_dot(&b, dir));
  vec3_add(&res, &res, &b);
  vec3_normalize(&res);
  return res;
}  

void add_circle(struct tangle_state *tangle,
		struct vec3d *center, struct vec3d *dir, double r,
		int Npoints)
{
  if(Npoints == 0)
      return;

  if(num_free_points(tangle) < Npoints)
    {
      expand_tangle(tangle, tangle->N + Npoints);
    }

  struct vec3d u, v, p, ptmp;
  int curr_point, last_point, first_point;
  struct vec3d zdir = perpendicular(dir);
  
  vec3_cross(&u, &zdir, dir);
  vec3_cross(&v, dir, &u);

  vec3_normalize(&u);
  vec3_normalize(&v);

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

int loop_injection = 0;
double loop_injection_frequency = 1;
void inject_loop(struct tangle_state *tangle, double t, double frequency)
{
  static double last_injection = 0;

  //ranges where to place the loop in the xy plane
  double XL, XH, YL, YH;

  //direction and center of the injected loop
  //inject in random direction pointing down-ish
  struct vec3d dir = vec3(2*(drand48() - 0.5), 2*(drand48()-0.5), -drand48());
  vec3_normalize(&dir);
  struct vec3d cent;

  //is it time yet to inject a loop?
  if(t - last_injection < 1/frequency)
      return;
  last_injection = t;

  XL = tangle->box.bottom_left_back.p[0];
  XH = tangle->box.top_right_front.p[0];

  YL = tangle->box.bottom_left_back.p[1];
  YH = tangle->box.top_right_front.p[1];

  //the injected loop is placed randomly in the top plane
  //of the domain box
  cent.p[2] = tangle->box.top_right_front.p[2];
  cent.p[0] = XL + (XH - XL)*drand48();
  cent.p[1] = YL + (YH - YL)*drand48();

  double D = fabs(XH - XL);
  double rmin = 0.05*D;
  double rmax = 0.25*D;

  double r = rmin + (rmax - rmin)*drand48();

  add_circle(tangle, &cent, &dir, r, 128);
}

int line_injection = 0;
int line_injection_n = 1;
double line_injection_frequency;
void inject_line_pairs(struct tangle_state *tangle, double t, double frequency)
{
  /*
   * Injects a pair of straignt vortex lines of opposite circulations oriented along z-axis.
   * Only makes sense for the 2-wall and periodic boundaries.
   */

  static double last_injection = 0;
  //is it time yet to inject the lines?
  if(t - last_injection < 1/frequency)
      return;
  last_injection = t;

  for(int k = 0; k < line_injection_n; ++k)
    {
      double XL = tangle->box.bottom_left_back.p[0];
      double XH = tangle->box.top_right_front.p[0];

      double YL = tangle->box.bottom_left_back.p[1];
      double YH = tangle->box.top_right_front.p[1];

      //get random positions
      double x1 = XL + drand48()*(XH-XL);
      double y1 = YL + drand48()*(YH-YL);
      double x2 = XL + drand48()*(XH-XL);
      double y2 = YL + drand48()*(YH-YL);

      add_line(tangle, x1, y1, +1, 20);
      add_line(tangle, x2, y2, -1, 20);
    }
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
  const double rmax = 0.25*D;

  for(int k = 0; k<N; ++k)
    {
      struct vec3d dir = vec3(drand48()-0.5, drand48()-0.5, drand48()-0.5);
      vec3_normalize(&dir);

      struct vec3d c;
      for(int j=0; j<3; ++j)
	c.p[j] = tangle->box.bottom_left_back.p[j] + drand48()*Ls[j];

      double r =  rmin + drand48()*(rmax - rmin);
      add_circle(tangle, &c, &dir, r, 256);
    }
}

void make_big_ring(struct tangle_state *tangle, double ring_r, int ring_N)
{
  /*
   * Creates a big vortex ring at [0,0,0] with thickness of 10% of radius
   */

  struct vec3d dir = vec3(0, 0, 1);


  for(int k=0; k<ring_N; ++k)
    {
      //get some randomness in [-1, 1]
      double sx = (drand48() - 0.5)*2;
      double sy = (drand48() - 0.5)*2;
      double sz = (drand48() - 0.5)*2;

      struct vec3d cent = vec3(sx,sy,sz);
      vec3_mul(&cent, &cent, 0.1*ring_r);

      double sr = (drand48() - 0.5)*2;
      double r = (1 + sr/10)*ring_r;

      add_circle(tangle, &cent, &dir, r, 64);
    }
}

void clip_at_wall(struct tangle_state *tangle)
{
  /*
   * Clips the tangle at the z-walls
   */
  //limits for unconstrained walls
  double llimit = -INFINITY;
  double ulimit = INFINITY;

  if(tangle->box.wall[Z_L] == WALL_MIRROR)
    llimit = tangle->box.bottom_left_back.p[2];
  if(tangle->box.wall[Z_H] == WALL_MIRROR)
    ulimit = tangle->box.top_right_front.p[2];

  enum POINT_STATES {
        OK=0,
        KILL,         //point to be clipped
        EDGE_FORWARD_L, //last point above the lower wall
        EDGE_REVERSE_L,  //first point above the lower wall
        EDGE_FORWARD_H, //last point below the upper wall
        EDGE_REVERSE_H  //first point below the upper wall
      };

  //recalculate is used as a tag for points to be clipped
  for(int k=0; k<tangle->N; ++k)
    tangle->recalculate[k] = OK;

  //first find all the points to be clipped
  for(int kk=0; kk<tangle->N; ++kk)
    {
      if(tangle->status[kk].status == EMPTY)
	continue;

      double zkk = tangle->vnodes[kk].p[2];
      if(zkk <= llimit || zkk >= ulimit)
	{
	  tangle->recalculate[kk] = KILL;
	  const int forward = tangle->connections[kk].forward;
	  const int reverse = tangle->connections[kk].reverse;

	  double zf = tangle->vnodes[forward].p[2];
	  double zr = tangle->vnodes[reverse].p[2];

	  if(zkk >= ulimit)
	    {
	      if(zf < ulimit)
		tangle->recalculate[forward] = EDGE_REVERSE_H;
	      if(zr < ulimit)
		tangle->recalculate[reverse] = EDGE_FORWARD_H;
	    }
	  if(zkk <= llimit)
	    {
	      if(zf > llimit)
		tangle->recalculate[forward] = EDGE_REVERSE_L;
	      if(zr > llimit)
		tangle->recalculate[reverse] = EDGE_FORWARD_L;
	    }
	}
    }

  //next handle edge points by projecting them on the wall
  //and pinning them on the lower or upper Z-wall
  for(int kk=0; kk < tangle->N; ++kk)
    {
      switch(tangle->recalculate[kk])
      {
	case EDGE_REVERSE_L:
	  tangle->connections[kk].reverse = -1;
	  tangle->vnodes[kk].p[2] = llimit;
	  tangle->status[kk].status = PINNED;
	  tangle->status[kk].pin_wall = Z_L;
	  break;
	case EDGE_FORWARD_L:
	  tangle->connections[kk].forward = -1;
	  tangle->vnodes[kk].p[2] = llimit;
	  tangle->status[kk].status = PINNED;
	  tangle->status[kk].pin_wall = Z_L;
	  break;
	case EDGE_REVERSE_H:
	  tangle->connections[kk].reverse = -1;
	  tangle->vnodes[kk].p[2] = ulimit;
	  tangle->status[kk].status = PINNED;
	  tangle->status[kk].pin_wall = Z_H;
	  break;
	case EDGE_FORWARD_H:
	  tangle->connections[kk].forward = -1;
	  tangle->vnodes[kk].p[2] = ulimit;
	  tangle->status[kk].status = PINNED;
	  tangle->status[kk].pin_wall = Z_H;
	  break;
	case KILL:
	  tangle->status[kk].status = EMPTY;
	  tangle->status[kk].pin_wall = NOT_A_FACE;
	  tangle->connections[kk].forward = -1;
	  tangle->connections[kk].reverse = -1;
	  break;
	default:
	  break;
      }
    }
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
  fprintf(stream, "\t%d\t%d\t%d\t%d\t%d", i,
	  tangle->connections[i].reverse,
	  tangle->connections[i].forward,
	  tangle->status[i].status,
	  tangle->status[i].pin_wall);
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
		  while(curr > 0 && tangle->connections[curr].reverse > 0)
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

  //stitch the last vortex together because the inner if never ran
  tangle->connections[pidx_start].reverse = pidx-1;
  tangle->connections[pidx-1].forward = pidx_start;

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

void add_line(struct tangle_state *tangle, double x, double y, int direction, int points)
{
  double zmin = tangle->box.bottom_left_back.p[2];
  double zmax = tangle->box.top_right_front.p[2];
  double dz = (zmax - zmin) / points;

  double zstart = direction > 0 ? zmin : zmax;
  double zend = direction > 0 ? zmax : zmin;

  struct vec3d s = vec3(x, y, zstart);
  struct vec3d sp = vec3(0, 0, direction); //tangent
  struct vec3d spp = vec3(0, 0, 0); //normal

  int new_pt = get_tangle_next_free(tangle);
  int last_pt;
  tangle->vnodes[new_pt] = s;
  tangle->tangents[new_pt] = sp;
  tangle->normals[new_pt] = spp;
  tangle->status[new_pt].status = PINNED;
  tangle->status[new_pt].pin_wall = direction > 0 ? Z_L : Z_H;
  tangle->connections[new_pt].reverse = -1;

  for(int k=1; k<points-1; ++k)
    {
      last_pt = new_pt;
      new_pt = get_tangle_next_free(tangle);
      s.p[2] = zstart + direction*k*dz;

      tangle->vnodes[new_pt] = s;
      tangle->tangents[new_pt] = sp;
      tangle->normals[new_pt] = spp;
      tangle->status[new_pt].status = FREE;
      tangle->status[new_pt].pin_wall = NOT_A_FACE;
      tangle->connections[new_pt].reverse = last_pt;
      tangle->connections[last_pt].forward = new_pt;
    }

  s.p[2] = zend;
  last_pt = new_pt;
  new_pt = get_tangle_next_free(tangle);
  tangle->vnodes[new_pt] = s;
  tangle->tangents[new_pt] = sp;
  tangle->normals[new_pt] = spp;
  tangle->status[new_pt].status = PINNED;
  tangle->status[new_pt].pin_wall = direction > 0 ? Z_H : Z_L;
  tangle->connections[new_pt].reverse = last_pt;
  tangle->connections[new_pt].forward = -1;
  tangle->connections[last_pt].forward = new_pt;
}

void random_straight_lines(struct tangle_state *tangle, int npairs, int points_per_line)
{
  //for the slab geometry, populate the computational box with straight vortices
  assert(tangle->box.wall[Z_L] == WALL_MIRROR && tangle->box.wall[Z_H] == WALL_MIRROR);

  double xmin = tangle->box.bottom_left_back.p[0];
  double xmax = tangle->box.top_right_front.p[0];
  double ymin = tangle->box.bottom_left_back.p[1];
  double ymax = tangle->box.top_right_front.p[1];

  for(int k=0; k < npairs; ++k)
    {
      double x1 = xmin + (xmax - xmin)*drand48();
      double y1 = ymin + (ymax - ymin)*drand48();
      double x2 = xmin + (xmax - xmin)*drand48();
      double y2 = ymin + (ymax - ymin)*drand48();
      add_line(tangle, x1, y1, +1, points_per_line);
      add_line(tangle, x2, y2, -1, points_per_line);
    }
}

double wall_dist(const struct tangle_state *tangle, int k, boundary_faces wall)
{
  /*
   * Checks whether a node k is closer than rdist to the wall.
   */
  int idx[6];
  idx[X_L] = idx[X_H] = 0;
  idx[Y_L] = idx[Y_H] = 1;
  idx[Z_L] = idx[Z_H] = 2;
  switch(wall)
  {
    case X_L:
    case Y_L:
    case Z_L:
      return tangle->vnodes[k].p[idx[wall]] - tangle->box.bottom_left_back.p[idx[wall]];

    case X_H:
    case Y_H:
    case Z_H:
      return tangle->box.top_right_front.p[idx[wall]] - tangle->vnodes[k].p[idx[wall]];

    default:
      error("wall_dist: unknown wall index %d\n", wall);
  }
  return -1;
}
