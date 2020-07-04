#include "debug.h"
#include "assert.h"
#include "brdf.h"
#include "litscene.h"
#include "polygon.h"
#include "sphere.h"
#include "simplecamera.h"
#include <GMmatlib.h>	
#include <math.h>
#include <float.h>


Ray SimpleCamera::StratifiedRandomRay(int i, int j, float k, float l, double offset)
{
	// [CS580 COMPUTER GRAPHICS] YOUR CODE HERE
	//==================================================
	//Skeleton for TASK2
	//Subtask 1
	//==================================================

	//(i,j) is the pixel location, 
	//and (k,l) is a stratified offset each for i and j. 
	//Finally offset is a coefficient for jittering.

	//Currently, we only sample a fixed point per a pixel. 
	//Implement the stratification and jittering technique in order to generate smooth shadow as shown in Figure 1.

	float kRand = k + frand(-offset, offset);
	float lRand = l + frand(-offset, offset);

	Point pixel_pos(Xmin + (i + kRand) * Width, Ymin + (j + lRand) * Height, 2.0);

	Ray ray;
	Point cop(0.0, 0.0, Zcop);

	ray.origin() = cop;
	ray.direction() = pixel_pos - cop;
	return ray;
}

Colour LitScene::renderPixel(int i, int j, SimpleCamera& TheCamera, int N_RAYS_PER_PIXEL)
{
	// [CS580 COMPUTER GRAPHICS] YOUR CODE HERE
	//==================================================
	//Skeleton for TASK2
	//Subtask 2
	//==================================================
	// (i,j): pixel location
	// You need to shoot N_RAYS_PER_PIXEL rays for a pixel (i,j) with the StratifiedRandomRay function, which you already implemented in SUBTASK(#1).
	// N_RAYS_PER_PIXEL = 64 - in this code
	Colour sum = Colour(0.0f, 0.0f, 0.0f);
	int cnt = 0;
	int subPixelWidth = sqrt(N_RAYS_PER_PIXEL);
	int subPixelHeight = N_RAYS_PER_PIXEL / subPixelWidth; // Case for not N_RAYS_PER_PIXEL being N^2
	for (int subi = 0; subi < subPixelWidth; subi++) {
		for (int subj = 0; subj < subPixelHeight; subj++) {
			float k = (0.5 + subi) / subPixelWidth;
			float l = (0.5 + subj) / subPixelHeight;
			double offset = 0.5 / max(subPixelHeight, subPixelWidth);
			sum = sum + tracePath(TheCamera.StratifiedRandomRay(i, j, k, l, offset));
		}
	}
	return sum / (subPixelHeight*subPixelWidth);
}


// ===========================================
// TASK 3-Area light sources
// ===========================================

float Polygon::TriangularSampling(Point& p, float s, float t, int Type = 0)
{
	// [CS580 COMPUTER GRAPHICS] YOUR CODE HERE
	//==================================================
	//Skeleton for TASK3
	//Subtask 1
	//==================================================
	// p is the sampled position of the light source. 
	// s and t are the coefficients for sampling. 
	// You can ignore the type parameter. 
	// Return value is the pdf for the p. 

	float x = (P[0].x() + P[1].x() + P[2].x()) / 3.0f;
	float y = (P[0].y() + P[1].y() + P[2].y()) / 3.0f;
	float z = (P[0].z() + P[1].z() + P[2].z()) / 3.0f;

	p = Point(x, y, z);
	float pdf = 1.0;
	return pdf;
}

float Polygon::RectangularSampling(Point& p, float s, float t, int Type = 0)
{
	// [CS580 COMPUTER GRAPHICS] YOUR CODE HERE
	//==================================================
	//Skeleton for TASK3
	//Subtask 2
	//==================================================
	// p is the sampled position of the light source. 
	// s and t are the coefficients for sampling. 
	// You can ignore the type parameter. 
	// Return value is the pdf for the p.

	float x = (P[0].x() + P[1].x() + P[2].x() + P[3].x()) / 4.0f;
	float y = (P[0].y() + P[1].y() + P[2].y() + P[3].y()) / 4.0f;
	float z = (P[0].z() + P[1].z() + P[2].z() + P[3].z()) / 4.0f;

	p = Point(x, y, z);
	float pdf = 1;
	return pdf;
}

bool Sphere::sample(Point& p, float& probability, const Point& from, float s, float t)
{
	// [CS580 COMPUTER GRAPHICS] YOUR CODE HERE
	//==================================================
	//Solution for TASK3
	//Subtask 3
	//==================================================
	// p is the sampled position of the light source. 
	// s and t are the coefficients for sampling. probability is the pdf for the p. 
	// from is the intersection point of a ray with an intersection object. 
	// *Reference: Section 3.2 ��Sampling Spherical Luminaries�� in ��Monte Carlo Techniques for Direct Lighting Calculations,�� ACM Transactions on Graphics, 1996

	probability = 1;
	p = (*this).Centre;

	return true;
}

// ===========================================
// TASK 4-Modified Phong BRDF-Diffuse, Specular
// ===========================================

Colour phongBRDF::brdf(
	const Point& x,				// surface point to evaluate brdf at
	const Vector& in,			// incoming vector
	const Vector& out,			// outgoing vector
	const Vector& Nx,			// normal at x (must be unit length)
	const float s
) const
{

	// [CS580 COMPUTER GRAPHICS] YOUR CODE HERE
	//==================================================
	//Skeleton for TASK4
	//Subtask 1
	//==================================================
	// A parameter s is the random variable in the range between 0 and 1. 
	// We use s in order to determine which term is employed from diffuse and specular BRDF.

	return m_kd / (float)PI;
}


// ===========================================
// TASK 5-Modified phong reflection
// ===========================================

Ray phongBRDF::reflection(const Ray& incoming,	// incoming ray
	const Vector& normal,	// surface normal
	const Point& point,	// intersection point with surface
	float s, float t,		// canonical random variables (on unit square)							  
	float& pdf
) const
{

	GPMatrix local2WC;
	GPVector3 U, N, V;
	GPVector	D;
	D.w = 1.0f;

	Ray r;
	r.origin() = point;

	float pd = m_kd.sum() / 3.0f;
	float ps = m_ks.sum() / 3.0f;
	float s2 = frand()*(pd + ps);
	if (s2 >= 0 && s2 < pd)
	{
		//==================================================
		//Diffuse part
		//==================================================

		D.x = cos(2.0f*(float)PI*s2)*sqrt(1 - t);
		D.y = sqrt(t);
		D.z = sin(2.0f*(float)PI*s2)*sqrt(1 - t);

		N.x = normal.x();
		N.y = normal.y();

		N.z = normal.z();

		// the vector is given in a right-handed orthonormal basis [U,normal,V], form that basis
		U = N.orthogonal();
		V = U % N;

		// now form a matrix that transforms the local coordinates to world coordinates
		local2WC.local2WCxForm(U, N, V);
		D = D * local2WC;
		r.direction().x() = D.x;
		r.direction().y() = D.y;
		r.direction().z() = D.z;

		pdf = (r.direction().normalised() ^ normal) / (float)PI;
		return r;
	}
	else if (s2 >= pd && s2 < pd + ps)
	{
		// [CS580 GLOBAL ILLUMINATION] YOUR CODE HERE
		//==================================================
		//TASK5
		//Specular part

	}

	pdf = 0;
	return r;
}

// ===========================================
// TASK 1-Direct illumiation
// TASK 6-Indirect illumination
// TASK 7-Russian roulette
// TASK 8-Refraction
// ===========================================

Colour LitScene::tracePath(const Ray& ray, GObject* object, Colour weightIn, int depth)
{
	float Lamda1 = frand(), Lamda2 = frand();   // random variables in cannonical 2d space in range of [0..1]

	Colour colourOut(0, 0, 0);

	// ----------------------------
	// check to see whether we should terminate here
	// ----------------------------
	if (depth >= m_maxBounces) return colourOut;
	if (weightIn.sum() < m_cutOffThreshold) return colourOut;

	Point ixp;		// intersection point
	Vector normal;	// normal at intersection point

	// ----------------------------
	// intersect ray with scene
	// ----------------------------
	GObject* hitObject = intersect(ray, ixp, normal);

	if (!hitObject) return colourOut;

	if (hitObject->isEmitter())
	{
		// ray struck an emitter - return emission.
		// you could possibly argue that there should be a slight chance that a ray is scattered 
		// off an emitter depending on its brdf - this is however usually ignored.
		float ndotl = normal ^ ray.direction().invert().normalised();
		if (ndotl > 0.0f)
			colourOut = hitObject->material().emission();
		return colourOut;
	}

	// [CS580 COMPUTER GRAPHICS] YOUR CODE HERE
	//==================================================
	//Skeleton for TASK1 
	//Subtask 1
	//==================================================
	// Direct illumination

	// If hitObject is normal object
	Colour colorAdd = Colour(0.0f, 0.0f, 0.0f);

	if (m_areaLightCount > 0) //If there is light source
	{
		//std::cout << "m_areaLightCount: " << m_areaLightCount << "\n";
		int theLight = int(frand() * m_areaLightCount);  // choose an emitter randomly

		Point lightPoint;
		float pdf;

		if (areaLightAt(theLight)->sample(lightPoint, pdf, ixp, Lamda1, Lamda2))
		{
			// generate ray in sampled direction
			Ray lightray;
			lightray.origin() = ixp;
			lightray.direction() = lightPoint - ixp;

			float ndotl = (normal ^ lightray.direction().normalised()); // dot product of normalised direction of light and surface
			if (ndotl > 0.0f) // If surface normal is same face as the light(if spotted positive face of surface)
			{
				Point intersectPoint;
				Vector normLight;
				if ((areaLightAt(theLight) == intersect(lightray, intersectPoint, normLight)))
				{
					float distance;
					float ndotlLight = lightray.direction().invert().normalised(distance) ^ normLight;
					if (ndotlLight > 0)
					{
						Colour brdf = hitObject->brdf()->brdf(ixp, ray.direction(), lightray.direction(), normal, hitObject->material().shininess());
						// Trace ray to emitter
						//std::cout << "areaLight: " << areaLightAt(theLight)->material().emission() << "\n";
						colorAdd = brdf * ndotl / pdf * (areaLightAt(theLight)->material().emission()*(ndotlLight / (distance*distance)));
					}
				}
			}
		}
	}

	// Indirect illumination
	// [CS580 GLOBAL ILLUMINATION] YOUR CODE HERE
	//==================================================
	//Task6, Task7, Task8
	//==================================================

	colourOut = colourOut + colorAdd;
	return colourOut;
}

// ===========================================
// Other functions
// ===========================================

Ray lambertianBRDF::reflection(
	const Ray& incoming,	// incoming ray
	const Vector& normal,	// surface normal
	const Point& point,	// intersection point with surface
	float s, float t,		// canonical random variables (on unit square)
	float& pdf
) const
{
	GPMatrix local2WC;
	GPVector3 U, N, V;
	GPVector	D;
	D.w = 1.0f;

	Ray r;
	r.origin() = point;

	D.x = cos(2.0f*(float)PI*s)*sqrt(1 - t);
	D.y = sqrt(t);
	D.z = sin(2.0f*(float)PI*s)*sqrt(1 - t);

	N.x = normal.x();
	N.y = normal.y();
	N.z = normal.z();

	// the vector is given in a right-handed orthonormal basis [U,normal,V], form that basis
	U = N.orthogonal();
	V = U % N;

	// now form a matrix that transforms the local coordinates to world coordinates
	local2WC.local2WCxForm(U, N, V);
	D = D * local2WC;
	r.direction().x() = D.x;
	r.direction().y() = D.y;
	r.direction().z() = D.z;

	pdf = (r.direction().normalised() ^ normal) / (float)PI;
	return r;
}

Colour lambertianBRDF::brdf(
	const Point& x,				// surface point to evaluate brdf at
	const Vector& in,			// incoming vector
	const Vector& out,			// outgoing vector
	const Vector& Nx,			// normal at x (must be unit length)
	const float s
) const
{
	return m_kd / (float)PI;
}