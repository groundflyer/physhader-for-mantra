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


// Geometry attenuation factor
//	nu - cosine of angle
float
gaf(float nu, alpha)
{
    return 2.0 / (1.0 + sqrt(1.0 + alpha * (1.0 / (nu*nu) - 1.0)));
}


// Distribution of visible normals
// slope formulation
vector2
get_slope(float theta, sx, _sy)
{
    float sy = _sy;

    if (theta < 0.0001)
	{
	    float r = sqrt(sx/(1.-sx));
	    float phi = 2. * PI * sy;
	    return set(r * cos(phi), r * sin(phi));
	}

    float tan_theta = tan(theta);
    float a = 1. / tan_theta;
    float G1 = 2. / (1. + sqrt(1. + 1./(a*a)));

    float A = 2. * sx/G1 - 1.;
    float A2 = A*A;
    float tmp = 1./(A2 - 1.);
    float B2 = tan_theta * tan_theta;
    float D = sqrt(B2 * tmp*tmp - (A2 - B2)*tmp);
    float slopex1 = tan_theta * tmp - D;
    float slopex2 = tan_theta * tmp + D;

    vector2 ret = 0;
    ret.x = (A < .0 || slopex2 > 1./tan_theta) ? slopex1 : slopex2;

    float S;
    if (sy > .5)
	{
	    S = 1.;
	    sy = 2. * (sy - .5);
	}
    else
	{
	    S = -1.;
	    sy = 2. * (.5 - sy);
	}

    float z = (sy*(sy*(sy * .27385 - .73369) + .46341)) /
	(sy*(sy*(sy * .093073 + .309420) - 1.) + .597999);
    ret.y = S * z * sqrt(1. + ret.x*ret.x);

    return ret;
}


// Distribution of visible normals microfacet
//	wi - incident direction
vector
microfacet(const vector _wi;
	   vector2 alpha;
	   float sx, sy)
{
    vector wi = _wi;
    wi.x *= alpha.x;
    wi.y *= alpha.y;

    wi = normalize(wi);

    float theta = .0;
    float phi = .0;

    if (wi.z < 0.99999)
	{
	    theta = acos(wi.z);
	    phi = atan2(wi.y, wi.x);
	}

    vector2 slope = get_slope(theta, sx, sy);

    float cosphi = cos(phi);
    float sinphi = sin(phi);
    slope *= set(cosphi, sinphi,
		 -sinphi, cosphi);

    slope *= alpha;

    return normalize(set(-slope.x, -slope.y, 1.));
}


// Get anisotropy rougnesses by anisotropy bias
vector2
anisorough(float alpha, bias)
{
    return set(alpha * (1. + bias), alpha * (1. - bias));
}


// incident projected roughness
float
anisorough_i(vector wi, tu, tv;
	     vector2 alpha)
{
    float cs = dot(wi, tu);
    float sn = dot(wi, tv);
    return sqrt(cs*cs*alpha.x*alpha.x + sn*sn*alpha.y);
}

// projected on surface rougness
float
anisorough_n(vector wm, wg, tu, tv;
	     float dotWmWg;
	     vector2 alpha)
{
    vector ph = normalize(wm - wg * dotWmWg);
    float
	cu = dot(ph, tu),
	cv = dot(ph, tv);
    return 1.0 / (cu*cu/alpha.x + cv*cv/alpha.y);
}


#endif	// __phy_microfacet__
