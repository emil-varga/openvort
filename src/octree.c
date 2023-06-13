/*
 * octree.c
 *
 *  Created on: Mar 26, 2019
 *      Author: emil
 */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "vortex_constants.h"
#include "octree.h"
#include "util.h"
#include "vec3_maths.h"

const int OCTREE_2D_CHILDREN_N = 4;

//create an empty octree with depth *depth* and room for *Ninit* vortex nodes
//at every tree node
struct octree* octree_init(int Ninit, int depth, int quadtree)
{
  struct octree *tree = (struct octree*)malloc(sizeof(struct octree));
  tree->quadtree = quadtree;

  const int children_N = tree->quadtree > 0 ? OCTREE_2D_CHILDREN_N : OCTREE_CHILDREN_N;

  tree->N=0;
  tree->N_total = Ninit;
  tree->node_ids = (int*)malloc(sizeof(int)*Ninit);
  if(depth > 1) {
    for(int k=0; k < children_N; ++k)
	  tree->children[k] = octree_init(Ninit, depth-1, quadtree);
  }
  else {
    for(int k=0; k < children_N; ++k)
      tree->children[k] = NULL;
  }
  return tree;
}

//free the entire structure whose root is *tree*
void octree_destroy(struct octree *tree)
{
  if(!tree)
    return;
  
  const int children_N = tree->quadtree > 0 ? OCTREE_2D_CHILDREN_N : OCTREE_CHILDREN_N;
  
  //we need to destroy it from the bottom up
  for(octree_child_idx child = 0; child < children_N; child++) {
    if(tree->children[child])
      octree_destroy(tree->children[child]);
  }
  free(tree);
}

//adds a node to the root of the tree
//this DOES NOT add the node recursively
struct octree *octree_add_node(struct octree *tree, int node_id)
{
  if(tree->N < tree->N_total) {
    tree->node_ids[tree->N++] = node_id;
    return tree;
  }

  tree->node_ids = (int*)realloc(tree->node_ids, 2*tree->N_total*sizeof(int));
  if(tree->node_ids) {
    tree->N_total *= 2;
    tree->node_ids[tree->N++] = node_id;
    return tree;
  }
  else
    return NULL;
}

int octree_inner_sort(struct octree *tree, const struct tangle_state *tangle);

//builds the octree from a tangle
struct octree *octree_build(const struct tangle_state *tangle, int quadtree)
{
  int total_children = 0;
  struct octree *tree = octree_init(tangle->N, 1, quadtree);
  //load all the non-empty points into the root of the tree
  for(int k=0; k < tangle->N; ++k) {
    if(tangle->status[k].status == EMPTY)
	    continue;
    if(!octree_add_node(tree, k))
	    error("failed to add node");
  }
  tree->tangle = tangle;
  tree->box = tangle->box;
  //sort the points in a recursive function into the children of the tree
  total_children += octree_inner_sort(tree, tangle);
  octree_update_means(tree, tangle);
  //printf("Total tree children: %d\n", total_children);
  return tree;
}

//sorting the nodes in the root into the children
int octree_inner_sort(struct octree *tree, const struct tangle_state *tangle)
{
  if(!tree || tree->N <= 1) //there is only one node, we do not need to split further
    return 0;
  
  if(min_box_size(&tree->box) < BH_grain*global_dl_max) {
    //the tree is already quite fine, don't split further
    //this is also necessary for consistency when using 2D tree
    return 0;
  }

  int new_children;
  if(tree->quadtree > 0)
    new_children = 4;
  else
    new_children = 8;
  octree_make_child_boxes(tree);
  const int children_N = tree->quadtree > 0 ? OCTREE_2D_CHILDREN_N : OCTREE_CHILDREN_N;
  for(int k=0; k < tree->N; ++k) {
    for(octree_child_idx child=0; child < children_N; child++) {
	    if(in_box(&tree->children[child]->box, &tangle->vnodes[tree->node_ids[k]])) {
        if(!octree_add_node(tree->children[child], tree->node_ids[k]))
		      error("failed to add node");
	      break;
	    }
	  }
  }

  for(octree_child_idx child=0; child < children_N; child++)
    new_children += octree_inner_sort(tree->children[child], tangle);
  
  //update the centers of mass, total circulation vectors and circulation tensors
  octree_update_means(tree, tangle);
  return new_children;
}

struct domain_box make_child_box(const struct domain_box *parent_box, octree_child_idx child_idx, int quadtree);
void octree_make_child_boxes(struct octree *tree)
{
  int Ninit = tree->N_total / 8;
  if(Ninit == 0) //less than 8 elements
    Ninit = 8;

  const int children_N = tree->quadtree > 0 ? OCTREE_2D_CHILDREN_N : OCTREE_CHILDREN_N;
  for(octree_child_idx child = 0; child < children_N; ++child) {
    struct domain_box child_box = make_child_box(&tree->box, child, tree->quadtree);
    tree->children[child] = octree_init(Ninit, 1, tree->quadtree); //do not create any depth here
    tree->children[child]->box = child_box;
    tree->children[child]->tangle = tree->tangle;
  }
}
//Helper for creating a child box. Gets the base dimensions from parent_box and returns the box
//corresponding to the child_idx.
struct domain_box make_child_box(const struct domain_box *parent_box, octree_child_idx child_idx, int quadtree)
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
      if(quadtree)
        error("no top boxes for quadtree");
      child_box.bottom_left_back = c;
      child_box.bottom_left_back.p[1] -= L.p[1]/2;
      child_box.top_right_front = parent_box->top_right_front;
      child_box.top_right_front.p[1] -= L.p[1]/2;
      break;
    case TLB: //top-left-back
      if(quadtree)
        error("no top boxes for quadtree");
      child_box.bottom_left_back = parent_box->bottom_left_back;
      child_box.bottom_left_back.p[2] += L.p[2]/2;
      child_box.top_right_front = c;
      child_box.top_right_front.p[2] += L.p[2]/2;
      break;
    case TRF: //top-right-front
      if(quadtree)
        error("no top boxes for quadtree");
      child_box.bottom_left_back = c;
      child_box.top_right_front = parent_box->top_right_front;
      break;
    case TRB: //top-right-back
      if(quadtree)
        error("no top boxes for quadtree");
      child_box.bottom_left_back = c;
      child_box.bottom_left_back.p[0] -= L.p[0]/2;
      child_box.top_right_front = parent_box->top_right_front;
      child_box.top_right_front.p[0] -= L.p[0]/2;
      break;
    default:
      error("make_child_box: incorrect child index");
      break;
  }

  if(quadtree)
    child_box.top_right_front.p[2] = parent_box->top_right_front.p[2];

  return child_box;
}

void octree_update_means(struct octree *tree, const struct tangle_state *tangle)
{
  if(!tree || tree->N == 0)
    return;

  tree->centre_of_mass = vec3(0,0,0);
  for(int k=0; k < tree->N; k++) {
    vec3_add(&tree->centre_of_mass, &tree->centre_of_mass, &tangle->vnodes[tree->node_ids[k]]);
  }
  vec3_mul(&tree->centre_of_mass, &tree->centre_of_mass, 1.0/tree->N);

  tree->total_circulation = vec3(0,0,0);
  tree->circ_tensor = mat3_null();
  for(int k=0; k < tree->N; k++) {
    int idx = tree->node_ids[k];
    int forward = tangle->connections[idx].forward;
    if(forward < 0)
      continue; //e.g., pinned points
    struct segment seg = seg_pwrap(&tangle->vnodes[idx], &tangle->vnodes[forward], &tangle->box);
    struct vec3d ds = segment_to_vec(&seg);
    vec3_add(&tree->total_circulation, &tree->total_circulation, &ds);

    struct vec3d sloc;
    vec3_sub(&sloc, &tangle->vnodes[idx], &tree->centre_of_mass);
    struct mat3d dM;
    vec3_outer(&dM, &ds, &sloc);
    mat3_add(&tree->circ_tensor, &tree->circ_tensor, &dM);
  }

  const int children_N = tree->quadtree ? OCTREE_2D_CHILDREN_N : OCTREE_CHILDREN_N;
  for(octree_child_idx child=0; child < children_N; child++)
    octree_update_means(tree->children[child], tangle);
}

//velocity calculation at point *r*
void octree_get_vs(const struct octree *tree, const struct vec3d *r, double resolution, struct vec3d *res, int skip)
{
  if(!tree || tree->N == 0) {
    //empty leaf, exit
    *res = vec3(0,0,0);
    return;
  }

  //check the resolution
  const double Lm = max_box_size(&tree->box);
  const double Lmin = min_box_size(&tree->box);
  struct segment seg = seg_pwrap(r, &tree->centre_of_mass, &tree->tangle->box);
  double d = segment_len(&seg);

  if(tree->N == 1 || Lmin < BH_grain) {
    //end leaf or the box is already small, integrate the BS directly
    *res = calculate_vs_shift(tree->tangle, *r, skip, NULL, tree->node_ids, tree->N);
    return;
  }

  if(d > 0 && Lm/d < resolution) {
    //calculate the velocity using the approximation
    struct vec3d v1;
    *res = vec3(0,0,0);
    struct vec3d R;
    double Rm;

    //first order
    vec3_sub(&R, r, &tree->centre_of_mass);
    Rm = vec3_d(&R);
    vec3_cross(&v1, &tree->total_circulation, &R);
    vec3_mul(&v1, &v1, 1/Rm/Rm/Rm);
    vec3_add(res, res, &v1);

    //second order, 1st term
    v1 = vec3(0,0,0);
    for(int i = 0; i<3; ++i) {
	    for(int j=i+1; j<3; ++j) {
        for(int k=j+1; k<3; ++k) {
		      v1.p[i] = epsilon[i][j][k]*tree->circ_tensor.m[k][j];
		    }
	    }
	  }
    vec3_mul(&v1, &v1, 1/Rm/Rm/Rm);
    vec3_add(res, res, &v1);

    // 2nd order, 2nd term
    mat3_vmul(&v1, &tree->circ_tensor, &R);
    vec3_cross(&v1, &v1, &R);
    vec3_mul(&v1, &v1, 1/Rm/Rm/Rm/Rm/Rm);
    vec3_add(res, res, &v1);

    vec3_mul(res, res, KAPPA/4/M_PI);

    return;
  }

  //resolution is not sufficient, we need to open the tree deeper
  struct vec3d partial_vs, total_vs;
  total_vs = vec3(0,0,0);
  const int children_N = tree->quadtree > 0 ? OCTREE_2D_CHILDREN_N : OCTREE_CHILDREN_N;
  for(int k = 0; k < children_N; ++k) {
      octree_get_vs(tree->children[k], r, resolution, &partial_vs, skip);
      vec3_add(&total_vs, &total_vs, &partial_vs);
  }
  *res = total_vs;
}
