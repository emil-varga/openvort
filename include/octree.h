/*
 * octree.h
 *
 *  Created on: Mar 26, 2019
 *      Author: emil
 */

#ifndef INCLUDE_OCTREE_H_
#define INCLUDE_OCTREE_H_

#include "vec3_maths.h"

typedef enum _octree_child_idx {
  BLF, BRF, BRB, BLB, //Bottom {Left/Right} {Front/Back}
  TLF, TRF, TRB, TLB  //Top {Left/Right} {Front/Back}
} octree_child_idx;

struct octree {
  int N_total;
  int N;
  struct vec3d *elements;
  struct vec3d centre_of_mass;
  struct vec3d total_circulation;
  struct domain_box box;
  struct octree *children[8];
};

struct octree* octree_init(int Ninit, int depth);
void octree_destroy(struct octree *tree);

int octree_add(struct octree *tree, const struct vec3d *v);
struct octree *octree_build(int N, const struct vec3d *vs);

/*helper functions*/
void octree_update_means(struct octree *tree);
void octree_create_children(struct octree *tree);

#endif /* INCLUDE_OCTREE_H_ */
