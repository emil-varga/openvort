/*
 * octree.c
 *
 *  Created on: Mar 26, 2019
 *      Author: emil
 */

#include <stdlib.h>
#include "octree.h"
#include "util.h"
#include "vec3_maths.h"

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

void octree_destroy(struct octree *tree)
{
  if(!tree)
    return;

  //we need to destroy it from the bottom up
  for(octree_child_idx child = 0; child < OCTREE_CHILDREN_N; child++)
    {
      if(tree->children[child])
	octree_destroy(tree->children[child]);
	free(tree->children[child]);
    }
  free(tree);
}

struct octree *octree_add_node(struct octree *tree, int node_id)
{
  if(tree->N < tree->N_total)
    {
      tree->node_ids[tree->N++] = node_id;
      return tree;
    }

  tree->node_ids = realloc(tree->node_ids, 2*tree->N_total);
  if(tree->node_ids)
    {
      tree->N_total *= 2;
      tree->node_ids[tree->N] = node_id;
      return tree;
    }
  else
    return NULL;
}

void octree_inner_sort(struct octree *tree, const struct tangle_state *tangle);
struct octree *octree_build(const struct tangle_state *tangle)
{
  struct octree *tree = octree_init(tangle->N, 1);
  //load up all the non-empty points
  for(int k=0; k < tangle->N; ++k)
    {
      if(tangle->status[k].status == EMPTY)
	continue;
      tree->node_ids[tree->N++] = k;
    }
  tree->box = tangle->box;
  //sort the points in a recursive function
  octree_inner_sort(tree, tangle);
  return tree;
}
void octree_inner_sort(struct octree *tree, const struct tangle_state *tangle)
{
  if(tree->N <= 1) //there is only one node, we do not need to split further
    return;

  octree_make_child_boxes(tree);
  for(int k=0; k < tree->N; ++k)
    {
      for(octree_child_idx child=0; child < OCTREE_CHILDREN_N; child++)
	{
	  if(in_box(&tree->box, &tangle->vnodes[tree->node_ids[k]]))
	    {
	      octree_add_node(tree->children[child], tree->node_ids[k]);
	      continue;
	    }
	}
    }
  for(octree_child_idx child=0; child < OCTREE_CHILDREN_N; child++)
    octree_inner_sort(tree->children[child], tangle);
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
  /*
   * Helper for creating a child box. Gets the base dimensions from parent_box and returns the box
   * corresponding to the child_idx.
   */
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
    case TLB: //top-left-back
      child_box.bottom_left_back = parent_box->bottom_left_back;
      child_box.bottom_left_back.p[2] += L.p[2]/2;
      child_box.top_right_front = c;
      child_box.top_right_front.p[2] += L.p[2]/2;
      break;
    case TRF: //top-right-front
      child_box.bottom_left_back = c;
      child_box.top_right_front = parent_box->top_right_front;
      break;
    case TRB: //top-right-back
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
