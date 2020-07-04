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

	Vector p0_p1 = P[1] - P[0];
	Vector p0_p2 = P[2] - P[0];

	// According to lectureNote15 page 30,
	// and page 18 of - Triangluar Sampling of  ¡°Monte Carlo Techniques for Direct Lighting Calculations, ¡± ACM Transactions on Graphics, 1996
	p = P[0] + p0_p1 * (t*sqrt(1 - s)) + p0_p2 * (1 - sqrt(1 - s));
	float pdf = 2.0 / (p0_p1*p0_p2).norm();
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

	Vector p0_p1 = P[1] - P[0];
	Vector p0_p3 = P[3] - P[0];

	p = P[0] + p0_p1 * s + p0_p3 * t;

	float pdf = 1.0 / (p0_p1*p0_p3).norm();
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
	// *Reference: Section 3.2 ¡®Sampling Spherical Luminaries¡¯ in ¡°Monte Carlo Techniques for Direct Lighting Calculations,¡± ACM Transactions on Graphics, 1996

	float distance = ((*this).Centre - from).norm();
	float tmp = ((*this).Radius / distance)*((*this).Radius / distance);
	double theta = acos(1 - s + t * sqrt(1 - tmp));
	double phi = 2 * (float)PI* t;
	Vector normHat = Vector(0, 1, 0);
	Vector wHat = normal(from).invert();
	Vector vHat = (wHat*normHat).normalised();
	Vector uHat = (vHat*wHat).normalised();
	double fx = cos(phi)*sin(theta);
	double fy = sin(phi)*sin(theta);
	double fz = cos(theta);

	Ray samplingRay;
	samplingRay.origin() = from;
	samplingRay.direction() = Vector(uHat.x()*fx + vHat.x()*fy + wHat.x()*fz, uHat.y()*fx + vHat.y()*fy + wHat.y()*fz, uHat.z()*fx + vHat.z()*fy + wHat.z()*fz);

	float factor;
	Colour color;
	(*this).intersect(samplingRay, factor, color);
	p = samplingRay.pointAt(factor);
	//std::cout << "samplingRay: " << samplingRay << "\n";
	//std::cout << "p: " << p << "\n";
	//std::cout << "From: " << from << "\n";
	//std::cout << "s and t : " << s << ", "<< t << "\n";

	Vector normSphere = normal(p);
	probability = (samplingRay.direction().invert() ^ normSphere) / (2 * (float)PI*factor*factor*(1 - sqrt(1 - tmp)));

	//probability = 1;
	//p = (*this).Centre;

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

	Colour diffuse = m_kd / (float)PI;

	Vector _in = in.normalised();
	Vector _out = out.normalised();
	Vector _Nx = Nx.normalised();
	Vector reflect = _Nx * (-2 * (_in ^ _Nx)) + _in;
	float cosAlpha = reflect.normalised() ^ _out;
	if (cosAlpha < 0.0f)
		cosAlpha = 0.0f;

	//specular term of modified phong brdf
	Colour specular = m_ks * ((s + 2) / (2.0*PI)) * pow(cosAlpha, s);

	return (diffuse + specular);
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
		//m_areaLightCount = 1 for all examples in this project
		int lightIdx = int(frand() * m_areaLightCount);  // select light source k

		Point lightSamplingPoint;
		float pdf;

		if (areaLightAt(lightIdx)->sample(lightSamplingPoint, pdf, ixp, Lamda1, Lamda2))
		{
			// ray in sampled direction
			Ray lightray;
			lightray.origin() = ixp;
			lightray.direction() = (lightSamplingPoint - ixp).normalised();

			float ndotl_surface = (normal ^ lightray.direction()); // dot product of normalised direction of light and surface
			if (ndotl_surface > 0.0f) // If surface normal is same face as the light(if spotted positive face of surface)
			{
				Point intersectPoint;
				Vector normLight;
				if ((areaLightAt(lightIdx) == intersect(lightray, intersectPoint, normLight)))
				{
					float distance = (lightSamplingPoint - ixp).norm();
					float ndotl_light = lightray.direction().invert() ^ normLight;
					if (ndotl_light > 0)
					{
						float s = hitObject->material().shininess();
						Colour brdf = hitObject->brdf()->brdf(ixp, ray.direction(), lightray.direction(), normal, s);
						// Trace ray to emitter
						Colour weightOut = (weightIn*brdf*ndotl_surface*(ndotl_light / (pdf*distance*distance)));
						//colorAdd = tracePath(lightray, hitObject, weightOut, depth + 1);
						//std::cout << "areaLight: " << areaLightAt(lightIdx)->material().emission() << "\n";
						if ((brdf.sum()*ndotl_surface*(1.0f / pdf)) > 3.0f) // If is over white(R+G+B>=3)
							colorAdd = (areaLightAt(lightIdx)->material().emission())*(ndotl_light / (distance*distance));
						else
							colorAdd = (areaLightAt(lightIdx)->material().emission()) * brdf * ndotl_surface *ndotl_light / (pdf*distance*distance);
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
