#ifndef TANGLE_H
#define TANGLE_H

#include <stdlib.h>
#include <stdint.h>

#include "vec3_maths.h"

typedef enum _ns_e {
  EMPTY = -1, //no point here
  FREE, //ordinary vortex point
  PINNED, //pinned on the wall
  PINNED_SLIP //pinned, but can slip
} node_status_t;

//status of the node and on which wall it is pinned
typedef struct _ns {
  node_status_t status;
  boundary_faces pin_wall;
} node_status;

#define PERIODIC_X 1
#define PERIODIC_Y 1<<1
#define PERIDOIC_Z 1<<2

struct neighbour_t {
  int forward;
  int reverse;
};

#define X_BSP_ID 0x00000000000FFFFF
#define Y_BSP_ID 0x000000FFFFF00000
#define Z_BSP_ID 0x0FFFFF0000000000

#define X_BSP_ID_OFFSET 0
#define Y_BSP_ID_OFFSET 20
#define Z_BSP_ID_OFFSET 40

/*
 * Structures describing the virtual tangles used for implementation of boundary conditions
 */

struct image_tangle {
  int shift[3]; //positive/negative number of shifts in the units of box size
  int reflect; //either a boundary_faces enum index or -1 for periodic
};

struct boundary_images {
  const struct image_tangle *images;
  int n;
};

//the actual image tangle configurations
//defined in boundary_images.c
extern const struct boundary_images open_boundaries;
extern const struct boundary_images periodic_6;
extern const struct boundary_images periodic_18;
extern const struct boundary_images periodic_26;

extern const struct boundary_images wall_1_6;
extern const struct boundary_images wall_1_18;
extern const struct boundary_images wall_1_26;


/*
 * The structure that holds all the tangle information
 */
struct tangle_state {
  struct vec3d *vnodes;
  struct vec3d *vnodes_new;
  struct vec3d *vels; //node velocities
  struct vec3d *vs; //superfluid velocity at the node
  struct vec3d *tangents;
  struct vec3d *normals;

  //flags that the properties of the points
  //need to be recalculated
  //currently used for not reconnecting a twice in a single pass
  int *recalculate;

  uint64_t *bsp_id; //for binary space partitioning -- id of the box

  struct neighbour_t *connections;

  node_status *status; //pinned, free, empty, etc.
  //the boundary conditions in the 6 cardinal directions
  struct domain_box box;
  struct boundary_images bimg;

  int N;
  int next_free;
  int total_free;
};

void create_tangle(struct tangle_state *tangle, size_t n);
void expand_tangle(struct tangle_state *tangle, size_t n);
void free_tangle(struct tangle_state *tangle);
struct vec3d step_node(const struct tangle_state *tangle, int i, int where);
int num_free_points(struct tangle_state *tangle);
int tangle_total_points(struct tangle_state *tangle);

//calculates the shifted r according to the image_tangle conf
struct vec3d shifted(const struct image_tangle *shift, const struct tangle_state *tangle,
		     const struct vec3d *r);
void enforce_boundaries(struct tangle_state *tangle);
void update_tangents_normals(struct tangle_state *tangle);
void update_velocities(struct tangle_state *tangle);
void update_velocity(struct tangle_state *tangle, int k);

/**
  @brief calculates the superfluid velocity

  Calculates and returns the superfluid velocity v_s induced by tangle at location r.
  Optionally disables one node given by skip. Use skip=-1 to not skip anything.
*/
struct vec3d calculate_vs(struct tangle_state *tangle, struct vec3d r, int skip);

void remesh(struct tangle_state *tangle, double min_dist, double max_dist);
void eliminate_small_loops(struct tangle_state *tangle, int loop_length);

int get_tangle_next_free(struct tangle_state *tangle);

static inline void update_tangle(struct tangle_state *tangle)
{
  update_tangents_normals(tangle);
  update_velocities(tangle);
}
#endif //TANGLE_H
