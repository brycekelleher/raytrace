#include <math.h>

bool SolveQuadratic(float *t0, float *t1, const float a, const float b, const float c)
{
	float d, rd, q;
	
	d = b * b - 4.0f * a * c;

	if(d < 0.0f)
		return false;

	rd = sqrtf(d);

	if(b < 0.0f)
		q = -0.5f * (b - rd);
	else
		q = -0.5f * (b + rd);

	*t0 = q / a;
	*t1 = c / q;

	return true;
}

// Theta is the sphere's elevation - using inclination gives a bad tangent
// vector. This means the sphere is parameterized as u [0..1] and v [-1 1]
// Clients of these functions will need to adjust their u, v values to match
// this parameterization.

void SpherePositionToUV(float *u, float *v, const float *xyz, float radius)
{
	float theta, phi;

	phi = atanf(xyz[1] / xyz[0]);
	theta = asinf(xyz[2] / radius);

	*u = phi / (2.0f * M_PI);
	*v = (2.0f * theta) / M_PI;
}

// Calculate a point on a unit sphere for a given uv
void SpherePosition(float *p, const float u, const float v)
{
	float theta, phi;

	phi = u * 2 * M_PI;
	theta = v * M_PI / 2;

	p[0] = cosf(phi) * cosf(theta);
	p[1] = sinf(phi) * cosf(theta);
	p[2] = sinf(theta);
}

// Calculate normalized tangents on a unit sphere for a given u, v
void SphereTangents(float *t0, float *t1, const float u, const float v)
{
	float theta, phi;

	phi = u * 2 * M_PI;
	theta = v * M_PI / 2;

	t0[0] = -sinf(phi) * cosf(theta);
	t0[1] = cosf(phi) * cosf(theta);
	t0[2] = 0.0f;

	t1[0] = -cosf(phi) * sinf(theta);
	t1[1] = -sinf(phi) * sinf(theta);
	t1[2] = cosf(theta);
}

bool SphereIntersect(const float *rayo, const float *rayd, const float radius, float *thit)
{
	// Compute intersection point
	// Compute u, v values
	// Compute Sphere Tangents
	float a, b, c, t0, t1;

	a = (rayd[0] * rayd[0]) + (rayd[1] * rayd[1]) + (rayd[2] * rayd[2]);
	b = 2.0f * ((rayd[0] * rayo[0]) + (rayd[1] * rayo[1]) + (rayd[2] * rayo[2]));
	c = ((rayo[0] * rayo[0]) + (rayo[1] * rayo[1]) + (rayo[2] * rayo[2])) - (radius * radius);

	if(!SolveQuadratic(&t0, &t1, a, b, c))
		return false;

	*thit = t0 < t1 ? t0 : t1;
	
	return true;
}

// Triangles have to be parametrized by hand (by assigning uv coordinates to the vertices)

// Triangle Intersection
// va = (v1 - v0) cross (v0 - v2) dot (v0 - a)
// v1 = ((v0 - o) cross d) dot (v0 - v2)
// v2 = ((v0 - o) cross d) dot (v1 - v0)
// v3 = ((v0 - o) cross d) dot (v2 - v1)
// u = v1 / (v1 + v2 + v3)
// v = v2 / (v1 + v2 + v3)
// w = v3 / (v1 + v2 + v3)

// fixme: Rename to ray triangle intersection test?
bool TriIntersect(const float *o, const float *d, const float *v0, const float *v1, const float *v2, float *uvw, float *t)
{
	float v0o[3], ov0[3], v1o[3], v2o[3], v1v0[3], v0v2[3];
	float v1o_cross_v2o[3], v2o_cross_v0o[3], v0o_cross_v1o[3], v1v0_cross_v0v2[3];
	float vola, vol0, vol1, vol2;
	bool hit;

	v0o[0] = v0[0] - o[0];
	v0o[1] = v0[1] - o[1];
	v0o[2] = v0[2] - o[2];
	ov0[0] = o[0] - v0[0];
	ov0[1] = o[1] - v0[1];
	ov0[2] = o[2] - v0[2];

	v1o[0] = v1[0] - o[0];
	v1o[1] = v1[1] - o[1];
	v1o[2] = v1[2] - o[2];

	v2o[0] = v2[0] - o[0];
	v2o[1] = v2[1] - o[1];
	v2o[2] = v2[2] - o[2];

	v1v0[0] = v1[0] - v0[0];
	v1v0[1] = v1[1] - v0[1];
	v1v0[2] = v1[2] - v0[2];

	v0v2[0] = v0[0] - v2[0];
	v0v2[1] = v0[1] - v2[1];
	v0v2[2] = v0[2] - v2[2];
	
	v1o_cross_v2o[0] = (v1o[1] * v2o[2]) - (v1o[2] * v2o[1]); 
	v1o_cross_v2o[1] = (v1o[2] * v2o[0]) - (v1o[0] * v2o[2]); 
	v1o_cross_v2o[2] = (v1o[0] * v2o[1]) - (v1o[1] * v2o[0]);
	
	v2o_cross_v0o[0] = (v2o[1] * v0o[2]) - (v2o[2] * v0o[1]); 
	v2o_cross_v0o[1] = (v2o[2] * v0o[0]) - (v2o[0] * v0o[2]); 
	v2o_cross_v0o[2] = (v2o[0] * v0o[1]) - (v2o[1] * v0o[0]);

	v0o_cross_v1o[0] = (v0o[1] * v1o[2]) - (v0o[2] * v1o[1]); 
	v0o_cross_v1o[1] = (v0o[2] * v1o[0]) - (v0o[0] * v1o[2]); 
	v0o_cross_v1o[2] = (v0o[0] * v1o[1]) - (v0o[1] * v1o[0]);

	vol0 = (v1o_cross_v2o[0] * d[0]) + (v1o_cross_v2o [1]* d[1]) + (v1o_cross_v2o[2] * d[2]);
	vol1 = (v2o_cross_v0o[0] * d[0]) + (v2o_cross_v0o [1]* d[1]) + (v2o_cross_v0o[2] * d[2]);
	vol2 = (v0o_cross_v1o[0] * d[0]) + (v0o_cross_v1o [1]* d[1]) + (v0o_cross_v1o[2] * d[2]);

	v1v0_cross_v0v2[0] = (v1v0[1] * v0v2[2]) - (v1v0[2] * v0v2[1]); 
	v1v0_cross_v0v2[1] = (v1v0[2] * v0v2[0]) - (v1v0[0] * v0v2[2]); 
	v1v0_cross_v0v2[2] = (v1v0[0] * v0v2[1]) - (v1v0[1] * v0v2[0]);

	vola = (v1v0_cross_v0v2[0] * ov0[0]) + (v1v0_cross_v0v2[1] * ov0[1]) + (v1v0_cross_v0v2[2] * ov0[2]);

	uvw[0] = vol0 / (vol0 + vol1 + vol2);
	uvw[1] = vol1 / (vol0 + vol1 + vol2);
	uvw[2] = vol2 / (vol0 + vol1 + vol2);
	*t = vola / (vol0 + vol1 + vol2);

	hit = (vol0 < 0 && vol1 < 0 && vol2 < 0) || (vol0 >= 0 && vol1 >= 0 && vol2 >= 0);

	return hit;
}

#if 0
void TriIntersect()
{
	// Test if collides (computes barycentrics)
	// Calculate intersection point from barycentrics
	// Calculate intersection normal
	// Calcuate uv values from barycentrics
	// Calcuate tangent space vectors
}
#endif

#include <stdlib.h>
#include <stdio.h>

void TestSphere()
{
	int i, j;

	for(i = 0; i <= 8; i++)
	{
		float u, v;
		float xyz[3], tangents[2][3];

		u = (float)i / 8.0f;
		v = (float)i / 8.0f;
		
		// convert from 0..1 to -1..1
		v = 2.0f * v - 1.0f;
		//v = -1.0f;

		SpherePosition(xyz, u, v);
		SphereTangents(tangents[0], tangents[1], u, v);

		fprintf(stdout, "%f %f %f %f %f %f\n", xyz[0], xyz[1], xyz[2], tangents[0][0], tangents[0][1], tangents[0][2]);
		fprintf(stdout, "%f %f %f %f %f %f\n", xyz[0], xyz[1], xyz[2], tangents[1][0], tangents[1][1], tangents[1][2]);
	}

}

void TestTriIntersect()
{
	float uvw[3], t;
	
	float vertices[] =
	{
		-1, -1, 1,
		 1, -1, 1,
		 0,  1, 1,
	};
	float rayo[] = { 0, 0, 0 };
	float rayd[] = { 1, 0, 0 };

	TriIntersect(rayo, rayd, &vertices[0], &vertices[3], &vertices[6], uvw, &t);
}

void TestTriIntersect2()
{
	float u, v, bary[3], t;
	extern float FoldedRadicalInverse(int n, int base);
	extern void CosineSampleHemisphere(float *vector, float u1, float u2);
	
	float rayo[3] = { 0, 0, 0 };
	float rayd[3];

#if 0
	float vertices[] =
	{
		-1, -1, 1,
		 1, -1, 1,
		 0,  1, 1,
	};
#endif
#if 1
	float vertices[] =
	{
		1, 0, 0,
		0, 1, 0,
		0, 0, 1,
	};
#endif

	for(int i = 0; i < 10000; i++)
	{
		u = FoldedRadicalInverse(i, 2);
		v = FoldedRadicalInverse(i, 3);
		
		CosineSampleHemisphere(rayd, u, v);

		if(TriIntersect(rayo, rayd, &vertices[0], &vertices[3], &vertices[6], bary, &t))
		{
			float hitpoint[3];

			hitpoint[0] = (bary[0] * vertices[0]) + (bary[1] * vertices[3]) + (bary[2] * vertices[6]);
			hitpoint[1] = (bary[0] * vertices[1]) + (bary[1] * vertices[4]) + (bary[2] * vertices[7]);
			hitpoint[2] = (bary[0] * vertices[2]) + (bary[1] * vertices[5]) + (bary[2] * vertices[8]);

			fprintf(stdout, "%f %f %f\n", hitpoint[0], hitpoint[1], hitpoint[2]);
		}
	}
}

void TestSphereIntersect()
{
	float u, v, thit;
	int i;
	extern float FoldedRadicalInverse(int n, int base);
	extern void CosineSampleHemisphere(float *vector, float u1, float u2);
	
	float rayo[3] = { 0, 0, -5.0f };
	float rayd[3];

	for(i = 0; i < 10000; i++)
	{
		u = FoldedRadicalInverse(i, 2);
		v = FoldedRadicalInverse(i, 3);
		
		CosineSampleHemisphere(rayd, u, v);

		if(SphereIntersect(rayo, rayd, 4.0f, &thit))
		{
			float hitpoint[3];

			hitpoint[0] = rayo[0] + (rayd[0] * thit);
			hitpoint[1] = rayo[1] + (rayd[1] * thit);
			hitpoint[2] = rayo[2] + (rayd[2] * thit);

			fprintf(stdout, "%f %f %f\n", hitpoint[0], hitpoint[1], hitpoint[2]);
		}
	}
}

#if 0
bool NewTriIntersect(const float *rayo, const float *rayd)
{
	// hit info structure
	float p[3];	// intersection point;
	float normal[3];
	float uv[2];	// uv parameter coordinates of hit point
	float dpdu[3];	// uv tangents
	float dpdv[3];

	if(TriIntersect(rayo, raydir, vertices[0], vertices[3], vertices[6], bary, &t))
	{
		// Check t value

		// Calculate the hit point
		p[0] = (bary[0] * vertices[0]) + (bary[1] * vertices[3]) + (bary[2] * vertices[6]);
		p[1] = (bary[0] * vertices[1]) + (bary[1] * vertices[4]) + (bary[2] * vertices[7]);
		p[2] = (bary[0] * vertices[2]) + (bary[1] * vertices[5]) + (bary[2] * vertices[8]);

		// Calculate the normal
		
		// Calculate the uv coordinates
		u[0] = (bary[0] * vertices[0]) + (bary[1] * vertices[3]) + (bary[2] * vertices[6]);
		v[1] = (bary[0] * vertices[1]) + (bary[1] * vertices[4]) + (bary[2] * vertices[7]);
		
		// Calculate the tangents
		dpdu[0] = (bary[0] * vertices[0]) + (bary[1] * vertices[3]) + (bary[2] * vertices[6]);
		dpdu[1] = (bary[0] * vertices[1]) + (bary[1] * vertices[4]) + (bary[2] * vertices[7]);
		dpdu[2] = (bary[0] * vertices[2]) + (bary[1] * vertices[5]) + (bary[2] * vertices[8]);
			
		dpdv[0] = (bary[0] * vertices[0]) + (bary[1] * vertices[3]) + (bary[2] * vertices[6]);
		dpdv[1] = (bary[0] * vertices[1]) + (bary[1] * vertices[4]) + (bary[2] * vertices[7]);
		dpdv[2] = (bary[0] * vertices[2]) + (bary[1] * vertices[5]) + (bary[2] * vertices[8]);
	}
}

bool NewSphereIntersect(const float *rayo, const float *rayd)
{
	float hitpoint[3];
	float tangents[2][3];
	float u, v;

	// Transform ray into sphere space
	
	SphereIntersect(rayo, raydir, hitpoint);

	// Calculate the uv coordinates
	SpherePositionToUV(hitpoint, &u, &v);

	// Calculate the normal
	
	// Calculate the tangents
	SphereTangents(tangents[0], tangents[1], u, v);
	
	// Scale the actual uv's into 0..1 space
	//v = 0.5f * v + 0.5f;
}
#endif

int main(int argc, char *argv[])
{
	//TestTriIntersect2();
	TestSphereIntersect();
}
