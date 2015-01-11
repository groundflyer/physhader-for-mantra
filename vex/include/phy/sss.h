// This may look like -*- C -*- code, but it is really Houdini Vex
// 
//	sss.h - Subsurface scattering routines.
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


#ifndef __phy_sss__
#define __phy_sss__


#include <pbr.h>
#include <phy/utils.h>
#include <phy/shading.h>


// Phase function
// Henyey-Greenstein distribution
//
// Parameters:
// 	* theta - cosine between incident and outgoing directions
// 	* g - average cosine of scattering
// Based on
// Henyey L. G. and Greenstein J. L.: Diffuse radiation in the galaxy.
float
phase(float theta, g)
{
    float g2 = g*g;
    return (1.0 - g2) / (1. + pow(1.0 + g2 - 2.0 * g * theta, 1.5));
}


// Fresnel angular moment
// 
// The approximation provided by
// d’Eon E., Irving G.: A quantized-diffusion model for
// rendering translucent materials.
void
fresnel_AM(float eta;
	   export float eval_2C1, eval_3C2)
{
    float
        eta2 = eta * eta,
        eta3 = eta2 * eta,
        eta4 = eta2 * eta2,
        eta5 = eta4 * eta;

    if (eta < 1.0)
	{
	    eval_2C1 = 0.919317 - 3.4793*eta + 6.75335*eta2
		- 7.80989*eta3 + 4.98554*eta4 - 1.36881*eta5;

	    eval_3C2 = 0.828421 - 2.62051*eta + 3.36231*eta2
		- 1.95284*eta3 + 0.236494*eta4 + 0.145787*eta5;
	}
    else
	{
	    eval_2C1 = -9.23372 + 22.2272*eta - 20.9292*eta2
		+ 10.2291*eta3 - 2.54396*eta4 + 0.254913*eta5;

	    eval_3C2 = -1641.1 + 135.926/eta3 - 656.175/eta2
		+ 1376.53/eta + 1213.67*eta - 568.556*eta2
		+ 164.798*eta3 - 27.0181*eta4 + 1.91826*eta5;
	}
}


// Irradiance for sss
vector
illum_surface(vector p, n;
	      float eta;
	      int sid;
	      int depth;
	      float depthimp)
{
    vector eval = .0;

    START_ILLUMINANCE;
    
    vector accum = .0;

    START_SAMPLING;
    SET_SAMPLE;

    float scale;
    vector lp, leval;
    mask = sample_light(lid, p, sample,
			Time, lp, leval, scale);

    cl = leval * scale / ALONE_VEC(leval);

    l = lp - p;
    cl *= shadow_light(lid, p, l, Time,
		       "noantialiasing", 1,
		       "N", n,
		       "SID", sid);

    if(mask & PBR_DIFFUSE_MASK)
	{
	    float kr, kt;
	    l = normalize(-l);
	    fresnel(l, n, eta, kr, kt);
	    accum += cl * kt;
	}

    END_LOOP; 			// SAMPLING

    eval += accum/samples;

    END_LOOP;			// ILLUMINANCE

    return eval;
}


// BSSRDF evaluator
// 
// Based on:
// 
// Eugene d'Eon: A Better Dipole (2012).
// 
// Ralf Habel, Per H. Christensen, Wojciech Jarosz:
// Classical and Improved Diffusion Theory for Subsurface Scattering (2013)
struct BSSRDF
{
    float g;			// Mean scattering cosine
    float eta;			// Relative IOR
    vector ca;			// Absorption cross-section
    vector cs;			// Scattering cross-section

    // Private
    vector muS_;		// Reduced scattering coefficient
    vector muT;			// Extinction coefficient
    vector muT_;		// Reduced extinction
    vector D;			// Diffusion coefficient
    vector muTR;		// Transport coefficient
    vector a_;			// Reduced scattering albedo
    vector Ce;			// Exitance parameter flux
    vector Cphi;		// Exitance parameter fluence

    float A_g;			// Reflection parameter
    float Zb;			// Imaginary depth

    float c1;			// Angular moments
    float c2;

    float mfp;			// Mean free path

    // Public

    void init(float _g, _eta;
	      vector _ca, _cs)
    {
	g = _g;
	eta = _eta;
	ca = _ca;
	cs = _cs;

	fresnel_AM(eta, c1, c2);

	muS_ = cs * (1.0 - g);
	muT = cs + ca;
	muT_ = muS_ + ca;

	A_g = (1.0 + c2) / (1.0 - c1);
	D = 0.333333333 * (1.0 / muT_ + ca / (muT_ * muT_));
	muTR = sqrt(ca / D);
	a_ = muS_ / muT_;
	Zb = 2.0 * max(D) * A_g;
	Cphi = 0.25 * (1.0 - c1);
	Ce = 0.5 * (1.0 - c2);

	mfp = 1. / max(muTR);
    }

    // BSSRDF approximation
    // from
    // Jensen H. W., Marschner S. R., Levoy M., Hanrahan P.:
    // A Practical Model for Subsurface Light Transport (2001)
    float approx()
    {
	float _a = max(a_);
	float ea = sqrt(3. * (1. - _a));

	return 0.5
	    * _a
	    * (1. + exp(-1.25 * A_g * ea))
	    * exp(-ea);
    }

    // Compute the BSSRDF
    // Xr - vector from incident point to real source
    // No - output normal
    vector eval(vector Xr, No)
    {
	vector Xrn = normalize(Xr);
	float kr, kt;
	vector nf = frontface(No, Xrn);

	float Zra = abs(dot(Xrn, nf));

	float
	    dR = length(Xr),
	    Zr = Zra * dR,
	    Zv = Zr + 2.0 * Zb;

	vector a_2 = a_ * a_;
	vector Xv = 2. * Zb * nf - Xr;

	float dV = length(Xv);
	vector
	    muTRdR = muTR * dR,
	    muTRdV = muTR * dV,
	    eSTDr = exp(-muTRdR),
	    eSTDv = exp(-muTRdV);

	vector Rphi =
	    Cphi
	    * a_2 / D
	    * (eSTDr / (dR + 1.)
	       - eSTDv / (dV + 1.));
	vector Re =
	    Ce
	    * a_2
	    * (Zr * eSTDr * (1. + muTRdR) / (dR*dR*dR + 1.)
	       + Zv * eSTDv * (1. + muTRdV) / (dV*dV*dV + 1.));

	fresnel(Xrn, nf, eta, kr, kt);

	vector _eval = (Rphi + Re) * kt;

	return _eval;
    }
}


// Compute real source position under point
vector
getRSP(vector p, n;
       float mfp)
{
    vector eval = p - n * mfp;

    // float bias;
    // renderstate("renderer:raybias", bias);
    // float depth = rayhittest(p, eval - p, bias);

    // if (depth > .0)
    // 	eval = p - n * depth * .99;

    return eval;
}


// Improved BSSRDF sampling strategy
// based on:
// King A., Kulla C., Conty A., Fajardo M.:
// BSSRDF Importance Sampling (2013)
vector
sampleSSS(float sx, sy, v)
{
    float theta = 2. * PI * sx;
    float z = sqrt(-2. * v * log(1. - v * (1. - exp(-0.5/v))));
    return set(sy * cos(theta), sy * sin(theta), z);
}

// Random tranformation matrix
matrix3
rand(float sx, sy)
{
    vector x = normalize(vector(rand(sx)) - 0.5);
    vector tmp = normalize(vector(rand(sy)) - 0.5);
    vector y = cross(x, tmp);
    vector z = cross(x, y);
    return set(x, y, z);
}


vector
raySSS(vector p, n;
       float eta, g;
       vector ca, cs;
       int samples, sid;
       string scope;
       int depth;
       float depthimp)
{
    BSSRDF bssrdf;
    bssrdf->init(g, eta, ca, cs);

    float v = max(bssrdf.muTR);
    float Rm = sqrt(12.46 / v);

    matrix3 basis;
    vector eval = .0;
    float pdf = .0;

    START_SAMPLING;

    vector pt = Rm * sampleSSS(sx, sy, v);
    vector dir = set(.0, .0, -1.);
    float dist = 2. * Rm * pt.z;

    basis = rand(sx, sy);

    pt = pt * basis + p;
    dir *= basis;

    vector hitP, hitN;

    if (trace(pt, dir, Time,
	      "samplefilter", "closest",
	      "scope", scope,
	      "maxdist", dist,
	      "P", hitP,
	      "N", hitN))
	{
	    float r = length(p - hitP);
	    hitN = normalize(hitN);
	    vector rsp = getRSP(hitP, hitN, bssrdf.mfp);
	    vector Xr = rsp - p;

	    float weight = exp(-r*r * .5 / bssrdf.mfp)
		/ (1. + abs(dot(n, hitN)));
	    vector hitC = illum_surface(hitP, hitN, eta, sid,
					depth, depthimp);
	    eval += weight * hitC * bssrdf->eval(Xr, n);
	    pdf += weight;
	}

    END_LOOP;

    return eval / pdf;
}


// Direct lighting for volume
vector
illum_volume(vector p, v;
	     float g;
	     int sid;
	     int depth;
	     float depthimp)
{
    vector eval = .0;

    START_ILLUMINANCE;
    START_SAMPLING;
    SET_SAMPLE;

    SAMPLE_LIGHT(p, v);

    if (g != 0)
	cl *= phase(dot(l,v), g);
    else
	cl *= .5;

    eval += cl;

    END_LOOP; 	// SAMPLING
    END_LOOP; 	// ILLUMINANCE

    return eval;
}


// Simple single scattering
vector
raymarch(vector p, v, ca, cs;
	 float g;
	 int samples, sid;
	 string scope;
	 int depth;
	 float depthimp)
{
    float
	pdf = .0,
	sc = (float)samples,
	maxdist = .0,
	raylength, bias;

    vector
	mu = cs + ca,
	eval = .0,
	pp, cl, l;

    float sigma = max(mu);
    vector _cs = cs / ALONE_VEC(cs);

    renderstate("renderer:raybias", bias);

    if (trace(p, v, Time,
	      "scope", scope,
	      "ray:length", raylength))
	maxdist = raylength;

    if (maxdist > bias * samples)
	{
	    START_SAMPLING;

	    if (maxdist > 1./sc)
		sx *= maxdist;

	    float weight = exp(-sigma * sx);

	    pdf += weight;

	    pp = p + v * sx;

	    cl = _cs * illum_volume(pp, v, g, sid,
				    depth, depthimp)
		* exp(-mu * sx);

	    eval += cl * weight;

	    END_LOOP;
	}

    if (pdf > .0)
    	eval /= pdf;

    return eval;
}


#endif
