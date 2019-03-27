/*
 * octree.c
 *
 *  Created on: Mar 26, 2019
 *      Author: emil
 */

#include <stdlib.h>
#include "octree.h"
#include "util.h"

struct octree* octree_init(int Ninit, int depth)
{
  struct octree *tree = (struct octree*)malloc(sizeof(struct octree));

  tree->N=0;
  tree->N_total = Ninit;
  tree->node_ids = (int*)malloc(sizeof(int)*Ninit);
  if(depth > 1)
    {
      for(int k=0; k < OCTREE_CHILDREN_N; ++k)
	tree->children[k] = octree_init(Ninit, depth-1);
    }
  return tree;
}

struct domain_box make_child_box(const struct domain_box *parent_box, octree_child_idx child_idx);
void octree_make_child_boxes(struct octree *tree)
{
  int Ninit = tree->N_total / 8;
  if(Ninit == 0) //less than 8 elements
    Ninit = 8;

  for(octree_child_idx child = 0; child < OCTREE_CHILDREN_N; ++child)
    {
      struct domain_box child_box = make_child_box(&tree->box, child);
      tree->children[child] = octree_init(Ninit, 1); //do not create any depth here
      tree->children[child]->box = child_box;
    }
}
struct domain_box make_child_box(const struct domain_box *parent_box, octree_child_idx child_idx)
{
  //box center and lengths
  struct vec3d c;
  struct vec3d L;

  vec3_add(&c, &parent_box->bottom_left_back, &parent_box->top_right_front);
  vec3_mul(&c, &c, 0.5);

  vec3_sub(&L, &parent_box->top_right_front, &parent_box->bottom_left_back);
  struct domain_box child_box;

  switch(child_idx)
  {
    //back->front -- increasing x
    //left->right -- increasing y
    //bottom->top -- increasing z
    case BLF: //bottom-left-front
      child_box.bottom_left_back = parent_box->bottom_left_back;
      child_box.bottom_left_back.p[0] += L.p[0]/2;
      child_box.top_right_front = c;
      child_box.top_right_front.p[0] += L.p[0]/2;
      break;
    case BLB: //bottom-left-back
      child_box.bottom_left_back = parent_box->bottom_left_back;
      child_box.top_right_front = c;
      break;
    case BRF: //bottom-right-front
      child_box.bottom_left_back = c;
      child_box.bottom_left_back.p[2] -= L.p[2]/2;
      child_box.top_right_front = parent_box->top_right_front;
      child_box.top_right_front.p[2] -= L.p[2]/2;
      break;
    case BRB: //bottom-right-back
      child_box.bottom_left_back = parent_box->bottom_left_back;
      child_box.bottom_left_back.p[1] += L.p[1]/2;
      child_box.top_right_front = c;
      child_box.top_right_front.p[1] += L.p[1]/2;
      break;
    case TLF: //top-left-front
      child_box.bottom_left_back = c;
      child_box.bottom_left_back.p[1] -= L.p[1]/2;
      child_box.top_right_front = parent_box->top_right_front;
      child_box.top_right_front.p[1] -= L.p[1]/2;
      break;
    case TLB:
      child_box.bottom_left_back = parent_box->bottom_left_back;
      child_box.bottom_left_back.p[2] += L.p[2]/2;
      child_box.top_right_front = c;
      child_box.top_right_front.p[2] += L.p[2]/2;
      break;
    case TRF:
      child_box.bottom_left_back = c;
      child_box.top_right_front = parent_box->top_right_front;
      break;
    case TRB:
      child_box.bottom_left_back = c;
      child_box.bottom_left_back.p[0] -= L.p[0]/2;
      child_box.top_right_front = parent_box->top_right_front;
      child_box.top_right_front.p[0] -= L.p[0]/2;
      break;
    default:
      error("make_child_box: incorrect child index");
      break;
  }

  return child_box;
}
