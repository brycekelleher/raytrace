#include "vec3.h"

typedef struct fsin_s
{
	vec2	framecoord;
	vec2	framesize;

} fsin_t;

typedef struct fsout_s;
{
	vec3	origin;
	vec3	dir;

} fsout_t;

static void frameshader()
{
	vec2 xy;
	xy = (in->framecoord / in->framesize);

	o.origin[0]  	= 0.0f;
	o.origin[1]	= 0.0f;
	o.origin[2]	= 0.0f;

	o.dir[0]	= xy[0] - 0.5f;
	o.dir[1]	= xy[1] - 0.5f;
	o.dir[2]	= 1.0f / (2.0f * tan(fov / 2.0f));
}

typedef struct rsin_s
{
	vec3	vertices[3];
	vec3	uvw;
	vec3	normal;
} rs_t;

static vec3 intersection()
{
	return (in->uvw[0] * in->vertices[0]) + (in->uvw[1] * in->vertices[1]) + (in->uvw[2] * in->vertices[2]);
}

static void rayshader()
{
	vec3 i	= intersection(in);
	vec3 n	= in->normal;
	vec3 l	= normalize(vec3(1, 1, 1));

	float ndotl = dot(n, l);
}
