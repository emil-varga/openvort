/*
 * octree.h
 *
 *  Created on: Mar 26, 2019
 *      Author: emil
 */

#ifndef INCLUDE_OCTREE_H_
#define INCLUDE_OCTREE_H_

#include "vec3_maths.h"
#include "tangle.h"

typedef enum _octree_child_idx {
  NA=-1,
  BLF, BLB, BRF, BRB, //Bottom {Left/Right} {Front/Back}
  TLF, TLB, TRF, TRB, //Top {Left/Right} {Front/Back}
  OCTREE_CHILDREN_N
} octree_child_idx;

extern const int OCTREE_2D_CHILDREN_N;

struct octree {
  int quadtree; //==0 for normal 3D octree, ==1 for 2D quadtree (do not split z-direction)
  int N_total; //total number of available spaces for vortex nodes in this tree node
  int N; //number of actual nodes in the this tree node
  int *node_ids;
  struct vec3d centre_of_mass;
  struct vec3d total_circulation;
  struct mat3d circ_tensor;
  struct domain_box box;
  struct octree *children[8];
  const struct tangle_state *tangle; //reference to the tangle used to build the tree
};

struct octree* octree_init(int Ninit, int depth, int quadtree);
void octree_destroy(struct octree *tree);

struct octree *octree_build(const struct tangle_state *tangle, int quadtree);

/*helper functions*/
octree_child_idx octree_find_child_index(const struct domain_box *box, const struct vec3d *r);
void octree_make_child_boxes(struct octree *tree);
void octree_update_means(struct octree *tree, const struct tangle_state *tangle);
void octree_create_children(struct octree *tree);

/*velocity calculation*/
void octree_get_vs(const struct octree *tree, const struct vec3d *r, double resolution, struct vec3d *res);

#endif /* INCLUDE_OCTREE_H_ */
