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

// derive the image count from the array itself, so the two cannot drift apart
#define IMAGES(arr) {.images = arr, .n = (int)(sizeof(arr) / sizeof((arr)[0]))}

const struct image_tangle periodic_z_open_xy_img[] = {{{0, 0, -1}, -1},
                                                      {{0, 0, 1}, -1}};

const struct image_tangle periodic_6_img[] = {
    {{-1, 0, 0}, -1}, {{0, -1, 0}, -1}, {{0, 0, -1}, -1},
    {{0, 0, 1}, -1},  {{0, 1, 0}, -1},  {{1, 0, 0}, -1}};

const struct image_tangle periodic_18_img[] = {
    {{-1, -1, 0}, -1}, {{-1, 0, -1}, -1}, {{-1, 0, 0}, -1}, {{-1, 0, 1}, -1},
    {{-1, 1, 0}, -1},  {{0, -1, -1}, -1}, {{0, -1, 0}, -1}, {{0, -1, 1}, -1},
    {{0, 0, -1}, -1},  {{0, 0, 1}, -1},   {{0, 1, -1}, -1}, {{0, 1, 0}, -1},
    {{0, 1, 1}, -1},   {{1, -1, 0}, -1},  {{1, 0, -1}, -1}, {{1, 0, 0}, -1},
    {{1, 0, 1}, -1},   {{1, 1, 0}, -1}};

const struct image_tangle periodic_26_img[] = {
    {{-1, -1, -1}, -1}, {{-1, -1, 0}, -1}, {{-1, -1, 1}, -1}, {{-1, 0, -1}, -1},
    {{-1, 0, 0}, -1},   {{-1, 0, 1}, -1},  {{-1, 1, -1}, -1}, {{-1, 1, 0}, -1},
    {{-1, 1, 1}, -1},   {{0, -1, -1}, -1}, {{0, -1, 0}, -1},  {{0, -1, 1}, -1},
    {{0, 0, -1}, -1},   {{0, 0, 1}, -1},   {{0, 1, -1}, -1},  {{0, 1, 0}, -1},
    {{0, 1, 1}, -1},    {{1, -1, -1}, -1}, {{1, -1, 0}, -1},  {{1, -1, 1}, -1},
    {{1, 0, -1}, -1},   {{1, 0, 0}, -1},   {{1, 0, 1}, -1},   {{1, 1, -1}, -1},
    {{1, 1, 0}, -1},    {{1, 1, 1}, -1}};

// wall in z, periodix in x, y
const struct image_tangle wall_2_4_img[] = {
    {{-1, 0, 0}, -1},                    // periodic in x
    {{1, 0, 0}, -1},  {{0, -1, 0}, -1},  // periodic in y
    {{0, 1, 0}, -1},  {{0, 0, -1}, Z_L}, // lower z-wall
    {{0, 0, 1}, Z_H}                     // upper z-wall
};

// as before, but also including the diagonal shifts
const struct image_tangle wall_2_6_img[] = {
    {{-1, 0, 0}, -1},                     // periodic in x
    {{1, 0, 0}, -1},   {{0, -1, 0}, -1},  // periodic in y
    {{0, 1, 0}, -1},   {{-1, -1, 0}, -1}, //-x-y
    {{-1, 1, 0}, -1},                     //-xy
    {{1, -1, 0}, -1},                     // x-y
    {{1, 1, 0}, -1},                      // xy
    {{0, 0, -1}, Z_L},                    // lower z-wall
    {{0, 0, 1}, Z_H}                      // upper z-wall
};

// as before, but also including the diagonal shifts also for the reflected
// tangles
const struct image_tangle wall_2_26_img[] = {
    {{-1, 0, 0}, -1},                       // periodic in x
    {{1, 0, 0}, -1},     {{0, -1, 0}, -1},  // periodic in y
    {{0, 1, 0}, -1},     {{-1, -1, 0}, -1}, //-x-y
    {{-1, 1, 0}, -1},                       //-xy
    {{1, -1, 0}, -1},                       // x-y
    {{1, 1, 0}, -1},                        // xy
    {{0, 0, -1}, Z_L},                      // lower z-wall
    {{1, 0, -1}, Z_L},                      // lower z-wall, x
    {{-1, 0, -1}, Z_L},                     // lower z-wall, -x
    {{0, 1, -1}, Z_L},                      // lower z-wall, y
    {{0, -1, -1}, Z_L},                     // lower z-wall, -y
    {{1, 1, -1}, Z_L},                      // lower z-wall, xy
    {{-1, 1, -1}, Z_L},                     // lower z-wall, -xy
    {{1, -1, -1}, Z_L},                     // lower z-wall, x-y
    {{-1, -1, -1}, Z_L},                    // lower z-wall, -x-y
    {{0, 0, 1}, Z_H},                       // upper z-wall
    {{1, 0, 1}, Z_H},                       // upper z-wall, x
    {{-1, 0, 1}, Z_H},                      // upper z-wall, -x
    {{0, 1, 1}, Z_H},                       // upper z-wall, y
    {{0, -1, 1}, Z_H},                      // upper z-wall, -y
    {{1, 1, 1}, Z_H},                       // upper z-wall, xy
    {{-1, 1, 1}, Z_H},                      // upper z-wall, -xy
    {{1, -1, 1}, Z_H},                      // upper z-wall, x-y
    {{-1, -1, 1}, Z_H},                     // upper z-wall, -x-y
};

const struct image_tangle wall_2_2_img[] = {
    {{-1, 0, 0}, -1},                     // periodic in x
    {{1, 0, 0}, -1},   {{0, -1, 0}, Y_L}, // front y wall
    {{0, 1, 0}, Y_H},                     // back y wall
    {{0, 0, -1}, Z_L},                    // lower z-wall
    {{0, 0, 1}, Z_H}                      // upper z-wall
};

const struct image_tangle wall_1_open_img[] = {{{0, 0, -1}, Z_L}};

const struct image_tangle wall_1_6_img[] = {{{-1, 0, 0}, -1},
                                            {{0, -1, 0}, -1},
                                            {{0, 0, -1}, Z_L},
                                            {{0, 0, 1}, -1},
                                            {{0, 1, 0}, -1},
                                            {{1, 0, 0}, -1}};

const struct image_tangle wall_1_18_img[] = {
    {{-1, -1, 0}, -1}, {{-1, 0, -1}, Z_L}, {{-1, 0, 0}, -1},  {{-1, 0, 1}, -1},
    {{-1, 1, 0}, -1},  {{0, -1, -1}, Z_L}, {{0, -1, 0}, -1},  {{0, -1, 1}, -1},
    {{0, 0, -1}, Z_L}, {{0, 0, 1}, -1},    {{0, 1, -1}, Z_L}, {{0, 1, 0}, -1},
    {{0, 1, 1}, -1},   {{1, -1, 0}, -1},   {{1, 0, -1}, Z_L}, {{1, 0, 0}, -1},
    {{1, 0, 1}, -1},   {{1, 1, 0}, -1},
};

const struct image_tangle wall_1_26_img[] = {
    {{-1, -1, -1}, Z_L}, {{-1, -1, 0}, -1}, {{-1, -1, 1}, -1},
    {{-1, 0, -1}, Z_L},  {{-1, 0, 0}, -1},  {{-1, 0, 1}, -1},
    {{-1, 1, -1}, Z_L},  {{-1, 1, 0}, -1},  {{-1, 1, 1}, -1},
    {{0, -1, -1}, Z_L},  {{0, -1, 0}, -1},  {{0, -1, 1}, -1},
    {{0, 0, -1}, Z_L},   {{0, 0, 1}, -1},   {{0, 1, -1}, Z_L},
    {{0, 1, 0}, -1},     {{0, 1, 1}, -1},   {{1, -1, -1}, Z_L},
    {{1, -1, 0}, -1},    {{1, -1, 1}, -1},  {{1, 0, -1}, Z_L},
    {{1, 0, 0}, -1},     {{1, 0, 1}, -1},   {{1, 1, -1}, Z_L},
    {{1, 1, 0}, -1},     {{1, 1, 1}, -1}};

const struct image_tangle channel_z_img[] = {
    {{-1, 0, 0}, X_L},                    // wall in x
    {{1, 0, 0}, X_H},  {{0, -1, 0}, Y_L}, // wall in y
    {{0, 1, 0}, Y_H},  {{0, 0, 1}, -1},   // periodic in z
    {{0, 0, -1}, -1}};

const struct boundary_images channel_z = IMAGES(channel_z_img);

const struct boundary_images open_boundaries = {.images = NULL, .n = 0};

const struct boundary_images periodic_z_open_xy = IMAGES(periodic_z_open_xy_img);

const struct boundary_images periodic_6 = IMAGES(periodic_6_img);

const struct boundary_images periodic_18 = IMAGES(periodic_18_img);

const struct boundary_images periodic_26 = IMAGES(periodic_26_img);

const struct boundary_images wall_1_open = IMAGES(wall_1_open_img);

const struct boundary_images wall_1_6 = IMAGES(wall_1_6_img);

const struct boundary_images wall_1_18 = IMAGES(wall_1_18_img);

const struct boundary_images wall_1_26 = IMAGES(wall_1_26_img);

const struct boundary_images wall_2_4 = IMAGES(wall_2_4_img);
const struct boundary_images wall_2_2 = IMAGES(wall_2_2_img);

const struct boundary_images wall_2_26 = IMAGES(wall_2_26_img);
