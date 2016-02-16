#include <stdio.h>
#include <stdlib.h>
#include "vec3.h"

float ScalarTriple(vec3 a, vec3 b, vec3 c)
{
	return Dot(a, Cross(b, c));
}

void IntersectTriangle(vec3 o, vec3 d, vec3 a, vec3 b, vec3 c)
{
	float	u, v, w;

	vec3 oa = a - o;
	vec3 ob = b - o;
	vec3 oc = c - o;

	u = ScalarTriple(d, oc, ob);
	v = ScalarTriple(d, oa, oc);
	w = ScalarTriple(d, ob, oa);

	float denom = 1.0f / (u + v + w);
	u *= denom;
	v *= denom;
	w *= denom;

	vec3 r = (u * a) + (v * b) + (w * c);
	printf("u=%f, v=%f, w=%f, r=(%f %f %f)\n", u, v, w, r[0], r[1], r[2]);
}

void TestIntersection(float *a, float *b, float *c)
{
	//vec3 va = vec3(a[0], a[1], a[2]);
	//vec3 vb = vec3(b[0], b[1], b[2]);
	//vec3 vc = vec3(c[0], c[1], c[2]);

	vec3 va = vec3(0, 0, 0);
	vec3 vb = vec3(2, 0, 0);
	vec3 vc = vec3(0, 1, 0);

	vec3 o = vec3(0.5, 0.5,  1);
	vec3 d = vec3(0, 0, -1);

	IntersectTriangle(o, d, va, vb, vc);
}

int main(int argc, char *argv[])
{
	TestIntersection(0, 0, 0);
}
