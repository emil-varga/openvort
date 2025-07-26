/*******************************************************************************
 * Copyright (C) 2019 Emil Varga <varga.emil@gmail.com>
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

#include "vortex_injection.h"
#include "vortex_utils.h"
#include "vortex_constants.h"
#include "math.h"

void inject_vortices(struct tangle_state *tangle, double t)
{
  if(global_loop_injection)
    inject_loop(tangle, t, global_loop_injection_frequency);
  if(global_line_injection)
    inject_line_pairs(tangle, t, global_line_injection_frequency, global_line_injection_n,
		      global_line_injection_polarized);
}

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

void inject_line_pairs(struct tangle_state *tangle, double t, double frequency,
		       int line_injection_n, int polarized)
{
  /*
   * Injects a pair of straight vortex lines of opposite circulations oriented along z-axis.
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

      //positions of the vortex pair
      double x1, y1, x2, y2;

      x1 = XL + drand48()*(XH-XL);
      x2 = XL + drand48()*(XH-XL);

      //get random positions
      if(polarized)
	{
	  y1 = YL + (1 + drand48())*(YH-YL)/2; //positive vortex in (0.5, 1)
	  y2 = YL + drand48()*(YH-YL)/2; //negative vortex in (0, 0.5)
	}
      else//unpolarized, fully random positions
	{
	  y1 = YL + drand48()*(YH-YL);
	  y2 = YL + drand48()*(YH-YL);
	}

      add_line(tangle, x1, y1, +1, 20);
      add_line(tangle, x2, y2, -1, 20);
    }
}
