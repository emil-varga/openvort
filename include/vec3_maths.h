#ifndef VEC3_MATHS_H
#define VEC3_MATHS_H

struct vec3d {
  double p[3];
};

struct segment {
  struct vec3d r1, r2;
};

typedef enum _wt {WALL_OPEN, WALL_PERIODIC, WALL_MIRROR} wall_type;

/*
 * Faces of the computational box. x low, x high etc.
 *
 * Functions dealing with boundaries in tangle.c depend on
 * the ordering of this enum.
 */
typedef enum _fc {
  X_L, X_H,
  Y_L, Y_H,
  Z_L, Z_H
} boundary_faces;

struct domain_box {
  struct vec3d bottom_left_back;
  struct vec3d top_right_front;
  wall_type wall[6];
};

extern const struct vec3d DIR_X, DIR_Y, DIR_Z;
extern const struct vec3d DIRS[3];

struct vec3d vec3(double x, double y, double z);
void vec3_assign(struct vec3d *v, double x, double y, double z);

struct domain_box make_box(struct vec3d bottom_left_front,
			   struct vec3d top_right_back,
			   wall_type wall[6]);

/*
 * open-space geometry
 */

double vec3_dot(const struct vec3d *u, const struct vec3d *v);
//normalized dot, cosine of the angle
double vec3_ndot(const struct vec3d *u, const struct vec3d *v);

void vec3_cross(struct vec3d *res,
		const struct vec3d *u, const struct vec3d *v);

void vec3_mul(struct vec3d *res,
	      const struct vec3d *u, double m);

void vec3_sub(struct vec3d *res,
	      const struct vec3d *u, const struct vec3d *v);

void vec3_add(struct vec3d *res,
	      const struct vec3d *u, const struct vec3d *v);

struct vec3d vec3_add2(const struct vec3d *u, const struct vec3d *v);

double vec3_d(const struct vec3d *u);

double vec3_dist(const struct vec3d *u, const struct vec3d *v);

void vec3_normalize(struct vec3d *v);

struct vec3d segment_to_vec(const struct segment *seg);
static inline double segment_len(const struct segment *seg)
{
  return vec3_dist(&seg->r1, &seg->r2);
}

/*
 * Periodic and mirror geometries
 */

struct segment seg_pwrap(const struct vec3d *r1, const struct vec3d *r2,
			 const struct domain_box *box);

struct vec3d mirror_shift(const struct vec3d *v, const struct domain_box *box,
			  boundary_faces wall);

struct vec3d periodic_shift(const struct vec3d *v, const struct domain_box *box,
			    boundary_faces wall);

#endif//VEC3_MATHS_H
