// This may look like -*- C -*- code, but it is really Houdini Vex
//
//	microfacet.h - BSDF routines for phySurface.
//	This is part of phyShader for Mantra.
//
// 
// Copyright (c) 2013-2015 Roman Saldygashev <sldg.roman@gmail.com>
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.


#ifndef __phy_microfacet__
#define __phy_microfacet__


#include <math.h>


// General parameters:
//	sigma - surface rougness
//	dotNH - dot product of normal and microfacet

// Next definitions are used in cvex shaders

#define GGG_REFL refl = 1

#define EVAL_GGG				\
    pdf = ct_ggg(sigma, dotNH, dotNV);		\
    eval = pdf;					\
    GGG_REFL

#define SAMPLE_GGG				\
    pdf = pdf_ggg(v, u, n, h, sigma);		\
    GGG_REFL

#define MAKE_MICROFACET				\
    vector h = microfacet(sigma, sx, sy);	\
    h *= set(tu, tv, n);

#define MAKE_ANISO_MICROFACET				\
    vector h = microfacet(sigmau, sigmav, sx, sy);	\
    h *= set(tu, tv, n);

#define ANISO_SIGMA_EVAL float sigma =			\
	sigmau == sigmav ? sigmau :			\
	anisorough(h, n, tu, tv, dotNH, sigmau, sigmav);
   
#define ANISO_SIGMA_PDF float sigma = avg(sigmau, sigmav);

#define REFLECTION_SAMPLE v = reflect(-u, h);

#define REFRACTION_EVAL				\
    vector tdir = -u;				\
    if (eta != 1.) tdir = refract(-u, n, eta);	\
    vector h = normalize(tdir + v);		\
    vector nb = frontface(n, u, n);		\
    float dotNH = dot(tdir, h);			\
    float dotNV = dot(nb, v);

#define REFRACTION_SAMPLE			\
    vector hn = n - h;				\
    hn = frontface(hn, u, hn);			\
    v = refract(normalize(hn-u), n, eta);


// Next routines based on:
// 
// Walter B., Marschner S. R., Li H., Torrance K. E.:
// Microfacet models for refraction through rough surfaces. (2007)
//
// Heitz E. and d’Eon E.:
// Importance Sampling Microfacet-Based BSDFs using
// the Distribution of Visible Normals. (2014)


// Modified GGX distribution
float
ggg(float dotNH, sigma)
{
    float D = sigma / (sigma + 1.0/(dotNH*dotNH) - 1.0);
    return D*D;
}


// Geometry attenuation factor
float
gaf(float nu, sigma)
{
    return 2.0 / (1.0 + sqrt(1.0 + sigma * (1.0 / (nu*nu) - 1.0)));
}


// Cook-Torrance (almost)
// Here is only shadowing GAF - masking term is in main routine
//	dotNL - dot product of normal and vector towards light
float
ct_ggg(float sigma, dotNH, dotNL)
{
    return ggg(dotNH, sigma) * gaf(dotNL, sigma);
}


// Sampling PDF
//	i - incoming light direction
//	o - outgoing light direction
//	n - normal
//	h - microfacet
float
pdf_ggg(vector i, o, n, h;
	float sigma)
{
    return abs(dot(i, h))
	* gaf(abs(dot(i, n)), sigma)
	* gaf(abs(dot(o, n)), sigma)
	/ (abs(dot(i, n)) + abs(dot(h, n)));
}


// Return sampled microfacet
// sx, sy - uniform random varieties on (0, 1)
// Isotropic case
vector
microfacet(float sigma, sx, sy)
{
    float tg = sigma * sqrt(sx) / sqrt(1. - sx);
    float theta = 2. * PI * sy;
    float x = tg * cos(theta);
    float y = tg * sin(theta);

    return normalize(set(x, y, 1.));
}

// Anisotropic case
vector
microfacet(float sigmau, sigmav, sx, sy)
{
    float
	_sx = sqrt(sx) / sqrt(1. - sx),
	theta = 2. * PI * sy;

    float x = sigmau * _sx * cos(theta);
    float y = sigmav * _sx * sin(theta);

    return normalize(set(x, y, 1.));
}


// Get anisotropy rougnesses by anisotropy bias
void
anisorough(float sigma, bias;
	   export float sigmau, sigmav)
{
    sigmau = sigma * (1. + bias);
    sigmav = sigma * (1. - bias);
}


// Compute rougness for current point
float
anisorough(vector ph, tu, tv;
	   float sigmau, sigmav)
{
    float
	cu = dot(ph, tu),
	cv = dot(ph, tv);
    return 1.0 / (cu*cu/sigmau + cv*cv/sigmav);
}

float
anisorough(vector h, n, tu, tv;
	   float dotNH, sigmau, sigmav)
{
    vector ph = normalize(h - n * dotNH);
    return anisorough(ph, tu, tv, sigmau, sigmav);
}


#endif	// __phy_microfacet__
