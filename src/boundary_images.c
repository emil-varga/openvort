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

#include "tangle.h"

/*
 * Constant boundary configuration declared in tangle.h
 */

const struct image_tangle periodic_6_img[] = {
    {{-1, 0, 0}, -1},
    {{0, -1, 0}, -1},
    {{0, 0, -1}, -1},
    {{0, 0, 1}, -1},
    {{0, 1, 0}, -1},
    {{1, 0, 0}, -1}
};

const struct image_tangle periodic_18_img[] = {
    {{-1, -1, 0}, -1},
    {{-1, 0, -1}, -1},
    {{-1, 0, 0}, -1},
    {{-1, 0, 1}, -1},
    {{-1, 1, 0}, -1},
    {{0, -1, -1}, -1},
    {{0, -1, 0}, -1},
    {{0, -1, 1}, -1},
    {{0, 0, -1}, -1},
    {{0, 0, 1}, -1},
    {{0, 1, -1}, -1},
    {{0, 1, 0}, -1},
    {{0, 1, 1}, -1},
    {{1, -1, 0}, -1},
    {{1, 0, -1}, -1},
    {{1, 0, 0}, -1},
    {{1, 0, 1}, -1},
    {{1, 1, 0}, -1}
};

const struct image_tangle periodic_26_img[] = {
    {{-1, -1, -1}, -1},
    {{-1, -1, 0}, -1},
    {{-1, -1, 1}, -1},
    {{-1, 0, -1}, -1},
    {{-1, 0, 0}, -1},
    {{-1, 0, 1}, -1},
    {{-1, 1, -1}, -1},
    {{-1, 1, 0}, -1},
    {{-1, 1, 1}, -1},
    {{0, -1, -1}, -1},
    {{0, -1, 0}, -1},
    {{0, -1, 1}, -1},
    {{0, 0, -1}, -1},
    {{0, 0, 1}, -1},
    {{0, 1, -1}, -1},
    {{0, 1, 0}, -1},
    {{0, 1, 1}, -1},
    {{1, -1, -1}, -1},
    {{1, -1, 0}, -1},
    {{1, -1, 1}, -1},
    {{1, 0, -1}, -1},
    {{1, 0, 0}, -1},
    {{1, 0, 1}, -1},
    {{1, 1, -1}, -1},
    {{1, 1, 0}, -1},
    {{1, 1, 1}, -1}
};

const struct image_tangle wall_2_4_img[] = {
    {{-1, 0, 0}, -1}, //periodic in x
    {{1, 0, 0}, -1},
    {{0, -1, 0}, -1}, //periodic in y
    {{0, 1, 0}, -1},
    {{0, 0, -1}, Z_L}, //lower z-wall
    {{0, 0, 1}, Z_H} //upper z-wall
};

const struct image_tangle wall_1_open_img[] = {
    {{0, 0, -1}, Z_L}
};

const struct image_tangle wall_1_6_img[] = {
    {{-1, 0, 0}, -1},
    {{0, -1, 0}, -1},
    {{0, 0, -1}, Z_L},
    {{0, 1, 0}, -1},
    {{1, 0, 0}, -1}
};

const struct image_tangle wall_1_18_img[] = {
    {{-1, -1, 0}, -1},
    {{-1, 0, -1}, Z_L},
    {{-1, 0, 0}, -1},
    {{-1, 0, 1}, -1},
    {{-1, 1, 0}, -1},
    {{0, -1, -1}, Z_L},
    {{0, -1, 0}, -1},
    {{0, -1, 1}, -1},
    {{0, 0, -1}, Z_L},
    {{0, 0, 1}, -1},
    {{0, 1, -1}, Z_L},
    {{0, 1, 0}, -1},
    {{0, 1, 1}, -1},
    {{1, -1, 0}, -1},
    {{1, 0, -1}, Z_L},
    {{1, 0, 0}, -1},
    {{1, 0, 1}, -1},
    {{1, 1, 0}, -1},
};

const struct image_tangle wall_1_26_img[] = {
    {{-1, -1, -1}, Z_L},
    {{-1, -1, 0}, -1},
    {{-1, -1, 1}, -1},
    {{-1, 0, -1}, Z_L},
    {{-1, 0, 0}, -1},
    {{-1, 0, 1}, -1},
    {{-1, 1, -1}, Z_L},
    {{-1, 1, 0}, -1},
    {{-1, 1, 1}, -1},
    {{0, -1, -1}, Z_L},
    {{0, -1, 0}, -1},
    {{0, -1, 1}, -1},
    {{0, 0, -1}, Z_L},
    {{0, 0, 1}, -1},
    {{0, 1, -1}, Z_L},
    {{0, 1, 0}, -1},
    {{0, 1, 1}, -1},
    {{1, -1, -1}, Z_L},
    {{1, -1, 0}, -1},
    {{1, -1, 1}, -1},
    {{1, 0, -1}, Z_L},
    {{1, 0, 0}, -1},
    {{1, 0, 1}, -1},
    {{1, 1, -1}, Z_L},
    {{1, 1, 0}, -1},
    {{1, 1, 1}, -1}
};

const struct boundary_images open_boundaries = {
    .images = NULL,
    .n = 0
};

const struct boundary_images periodic_6 = {
   .images = periodic_6_img,
   .n = 6
};

const struct boundary_images periodic_18 = {
   .images = periodic_18_img,
   .n = 18
};

const struct boundary_images periodic_26 = {
   .images = periodic_26_img,
   .n = 26
};

const struct boundary_images wall_1_open = {
   .images = wall_1_open_img,
   .n = 1
};

const struct boundary_images wall_1_6 = {
   .images = wall_1_6_img,
   .n = 6
};

const struct boundary_images wall_1_18 = {
   .images = wall_1_18_img,
   .n = 18
};

const struct boundary_images wall_1_26 = {
   .images = wall_1_26_img,
   .n = 26
};

const struct boundary_images wall_2_4 = {
    .images = wall_2_4_img,
    .n = 6
};
