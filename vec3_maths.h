#ifndef VEC3_MATHS_H
#define VEC3_MATHS_H

struct vec3d {
  double p[3];
};

struct vec3d vec3(double x, double y, double z);

void vec3_assign(struct vec3d *v, double x, double y, double z);

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

double vec3_d(const struct vec3d *u);

double vec3_dist(const struct vec3d *u, const struct vec3d *v);

#endif//VEC3_MATHS_H
