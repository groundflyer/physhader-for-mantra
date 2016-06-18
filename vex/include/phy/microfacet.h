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
//	alpha - surface rougness
//	dotWmWg - dot product of normal and microfacet

// Next definitions are used in cvex shaders


// GGX distribution
float
ggg(float dotWmWg, alpha)
{
    float D = alpha / (alpha + 1.0/(dotWmWg*dotWmWg) - 1.0);
    return D*D;
}

float
ggxalbedo(float alpha)
{
    float tmp = alpha + 1.;
    return 0.25 * (3.0 * alpha) / (tmp*tmp);
}

// Geometry attenuation factor
//	nu - cosine of angle
float
gaf(float nu, alpha)
{
    return 2.0 / (1.0 + sqrt(1.0 + alpha * (1.0 / (nu*nu) - 1.0)));
}


// Here is only shadowing GAF - masking term is in main routine
//	dotWgWi - dot product of normal and vector towards light
float
ct_ggg(float alpha, dotWmWg, dotWgWi)
{
    return ggg(dotWmWg, alpha) * gaf(dotWgWi, alpha);
}


// Sampling PDF
//	wi - incoming light direction
//	wo - outgoing light direction
//	wg - normal
//	wm - microfacet
float
pdf_ggg(vector wi, wo, wg, wm;
	float alpha)
{
    return abs(dot(wi, wm))
	* gaf(abs(dot(wi, wg)), alpha)
	* gaf(abs(dot(wo, wg)), alpha)
	/ (abs(dot(wi, wg)) + abs(dot(wm, wg)));
}


// Return sampled microfacet
// sx, sy - uniform random varieties on (0, 1)
// Isotropic case
vector
microfacet(float alpha, sx, sy)
{
    float tg = alpha * sqrt(sx) / sqrt(1. - sx);
    float theta = 2. * PI * sy;
    float x = tg * cos(theta);
    float y = tg * sin(theta);

    return normalize(set(x, y, 1.));
}

// Anisotropic case
vector
microfacet(float alphau, alphav, sx, sy)
{
    float
	_sx = sqrt(sx) / sqrt(1. - sx),
	theta = 2. * PI * sy;

    float x = alphau * _sx * cos(theta);
    float y = alphav * _sx * sin(theta);

    return normalize(set(x, y, 1.));
}


// Get anisotropy rougnesses by anisotropy bias
void
anisorough(float alpha, bias;
	   export float alphau, alphav)
{
    alphau = alpha * (1. + bias);
    alphav = alpha * (1. - bias);
}


// Compute rougness for current point
float
anisorough(vector ph, tu, tv;
	   float alphau, alphav)
{
    float
	cu = dot(ph, tu),
	cv = dot(ph, tv);
    return 1.0 / (cu*cu/alphau + cv*cv/alphav);
}

float
anisorough(vector wm, wg, tu, tv;
	   float dotWmWg, alphau, alphav)
{
    vector ph = normalize(wm - wg * dotWmWg);
    return anisorough(ph, tu, tv, alphau, alphav);
}


#endif	// __phy_microfacet__
