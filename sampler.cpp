#include <stdio.h>
#include <stdlib.h>
#include <math.h>

float RandomFloat()
{
	return (float)rand() / RAND_MAX;
}

float RandomFloat48()
{
	return (float)drand48();
}

// a 2d halton sequence would be x = ((RadicalInverse(2), RadicalInverse(3))
// a 2d hammersley sequence would be x = ((i / N), RadicalInverse(2))
// Both use prime numbers for the base as dimensionality increases

float RadicalInverse(int n, int base)
{
	float val, invbase, invbi;

	val = 0;
	invbase = 1.0f / base;
	invbi = invbase;

	while(n > 0)
	{
		val += (n % base) * invbi;
		n /= base;
		invbi *= invbase;
	}

	return val;
}

float FoldedRadicalInverse(int n, int base)
{
	float val, invbase, invbi;
	int modoffset, digit;

	val = 0;
	invbase = 1.0f / base;
	invbi = invbase;
	modoffset = 0;

	while(val + base * invbi != val)
	{
		digit = ((n + modoffset) % base);
		val += digit * invbi;
		n /= base;
		invbi *= invbase;
		modoffset++;
	}

	return val;
}

void UniformSampleHemisphere(float *vector, float u1, float u2)
{
	float r, phi;
	
	r = sqrtf(1.0f - u1 * u1);
	phi = 2.0f * M_PI * u2;
	vector[0] = r * cosf(phi);
	vector[1] = r * sinf(phi);
	vector[2] = u1;
}

void CosineSampleHemisphere(float *vector, float u1, float u2)
{
	float r, theta;

	r = sqrtf(u1);
	theta = 2 * M_PI * u2;
	vector[0] = r * cosf(theta);
	vector[1] = r * sinf(theta);
	vector[2] = sqrtf(1.0f - u1);
}

bool PointInTriangle2D(float *p, float *xy0, float *xy1, float *xy2)
{
	float e0[2], e1[2], e2[2], area[3];

	e0[0] = xy0[0] - p[0];
	e0[1] = xy0[1] - p[1];
	e1[0] = xy1[0] - p[0];
	e1[1] = xy1[1] - p[1];
	e2[0] = xy2[0] - p[0];
	e2[1] = xy2[1] - p[1];

	area[0] = (e0[0] * e1[1]) - (e1[0] * e0[1]);
	area[1] = (e1[0] * e2[1]) - (e2[0] * e1[1]);
	area[2] = (e2[0] * e0[1]) - (e0[0] * e2[1]);

	return (area[0] > 0.0f) && (area[1] > 0.0f) && (area[2] > 0.0f);
}

void TestPointInTriangle()
{
	float xy0[2] = { 0, 0 };
	float xy1[2] = { 20, 5 };
	float xy2[2] = { 0, 10 };
	float p[2] = { 5, 5 };

	PointInTriangle2D(p, xy0, xy1, xy2);
}

#if 0
int main(int argc, char *argv[])
{
	int i;
	FILE *fp;
	float x, y;

	int numsamples = 256;

	for(i = 0; i < numsamples; i++)
	{
		float x, y;
		//fprintf(stdout, "%f\n", RadicalInverse(i, 2));
		
		x = FoldedRadicalInverse(i, 2);
		y = FoldedRadicalInverse(i, 3);
		x = RandomFloat48();
		y = RandomFloat48();
		x = (float)i / numsamples;
		y = RadicalInverse(i, 2);
		x = FoldedRadicalInverse(i, 2);
		y = FoldedRadicalInverse(i, 3);
		//fprintf(stdout, "%f %f\n", x, y);
		//fprintf(stdout, "%f %f\n", ((float)i / (float)8), x);
		
		float vector[3];
		UniformSampleHemisphere(vector, x, y);
		fprintf(stdout, "%f %f %f\n", vector[0], vector[1], vector[2]);
	}

	TestPointInTriangle();

	return 0;
}
#endif
