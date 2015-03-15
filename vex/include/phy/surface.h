// This may look like -*- C -*- code, but it is really Houdini Vex
//
//	surface.h - phySurface for Mantra,
//	Physical based, easy to use, compact surface ubershader.
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


#ifndef __physurface__
#define __physurface__


#include <pbr.h>
#include <phy/utils.h>
#include <phy/spectrum.h>
#include <phy/microfacet.h>
#include <phy/volume.h>


// Treat surface smooth if roughness below this value
#define SMOOTH_THRESHOLD	.00005

#define MAX_ROUGH	0.3


// Direct lighting
void
illum_surface(vector p, pTRN, pSSS;
	      vector nfN, nbN;
	      vector v;
	      int thick;
	      float thickness;
	      float depthimp;
	      int sid;
	      int depth;
	      int enableDFS, enableSPC, enableTRN, enableSSS;
	      bsdf f_DFS, f_SPC, f_TRN, f_SSS;
	      vector factorDFS, factorSPC, factorTRN, factorSSS;
	      export vector dfs, spc, trn, sss)
{
    vector
	_dfs = .0,
	_spc = .0,
	_trn = .0,
	_sss = .0;

    START_ILLUMINANCE;

    vector _tmpDFS = .0;
    vector _tmpSPC = .0;
    vector _tmpTRN = .0;
    vector _tmpSSS = .0;

    START_SAMPLING("nextpixel");
    SET_SAMPLE;

    SAMPLE_LIGHT(p, nfN);

    // Diffuse reflection
    if (enableDFS && (mask & PBR_DIFFUSE_MASK))
	_tmpDFS += cl * EVAL_BSDF(f_DFS);

    // Specular reflection
    if (enableSPC && (mask & PBR_REFLECT_MASK))
	_tmpSPC += cl * EVAL_BSDF(f_SPC);

    // Refraction
    if (enableTRN)
	{
	    if (thick)
		SAMPLE_LIGHT(pTRN, nbN);

	    if (mask & PBR_REFRACT_MASK)
		_tmpTRN += cl * EVAL_BSDF(f_TRN);
	}

    // SSS of thin objects
    if (enableSSS)
	{
	    if (thick)
		SAMPLE_LIGHT(pSSS, nbN);

	    if (mask & PBR_DIFFUSE_MASK)
		_tmpSSS += cl * EVAL_BSDF(f_SSS);
	}

    END_LOOP;			// SAMPLING

    dfs = factorDFS * _tmpDFS / samples;
    spc = factorSPC * _tmpSPC / samples;
    trn = factorTRN * _tmpTRN / samples;
    sss = factorSSS * _tmpSSS / samples;

    if (!depth)
	{
	    vector all_comp[] = array(dfs, spc, trn, sss);
	    storelightexport(getlightname(lid), "all_comp", all_comp);
	}

    _dfs += dfs;
    _spc += spc;
    _trn += trn;
    _sss += sss;

    END_LOOP;			// ILLUMINANCE

    dfs = _dfs;
    spc = _spc;
    trn = _trn;
    sss = _sss;
}


// Complex Fresnel term
// n = etai/etat
// k - extinction coefficient
float
cfresnel(vector v, normal;
	 float n, k)
{
    float u = dot(v, normal);
    float n2 = n * n;
    float k2 = k * k;
    float u2 = u * u;
    float nk4 = 4. * n2 * k2;
    float n2minusk2 = n2 - k2;
    float u2minus1 = u2 - 1.;
    float tmp01 = n2minusk2 + u2minus1;
    float tmp2 = tmp01 * tmp01;
    float tmp2nk4 = tmp2 + nk4;
    float tmp3 = tmp2nk4 + n2minusk2 + u2minus1;
    float a2 = sqrt(tmp3) / 2.;
    float tmp4 = tmp2nk4 - n2 + k2 - u2 + 1.;
    float b2 = sqrt(tmp4) / 2.;
    float a = sqrt(a2);
    float aminusu = a - u;
    float aplusu = a + u;
    float F1 = (aminusu*aminusu + b2) / (aplusu*aplusu + b2) / 2.;
    float u1 = 1. / u;
    float aminusu1 = aminusu + u1;
    float aplusu1 = aplusu - u1;
    float F2 = (aplusu1*aplusu1 + b2) / (aminusu1*aminusu1 + b2) + 1.;
    return F1 * F2;
}


// Trace reflection/refraction
// Style choices:
//	1 - Gather
//	2 - Occlusion
#define TRACE_FLAGS(N) "scope", scope,			\
	"samplefilter", sfilter,			\
	"maxdist", maxdist,				\
	"samples", N,					\
	"raystyle", raystyle

#define TRACE_GATHER trace(p, dir, Time,		\
			   TRACE_FLAGS(1),		\
			   variable, hitCf,		\
			   "variancevar", variable)

#define TRACE_OCCL trace(p, dir, Time, TRACE_FLAGS(1))

#define TRACE_ABSRP trace(p, dir, Time,			\
			  TRACE_FLAGS(1),		\
			  variable, hitCf,		\
			  "ray:length", raylength)

#define INIT_TRACE vector hitCf = .0;			\
    vector eval = .0;					\
    int mask = PBR_REFRACT_MASK;			\
    if (raystyle == "reflect") mask = PBR_REFLECT_MASK; \
    string sfilter = oblend ? "opacity" : "closest";	\
    int doabs = (max(absty) > .0);			\
    float raylength = .0;

// Delta function case
vector
raytrace(vector p, dir;
	 float maxdist;
	 int oblend, style;
	 string raystyle, scope, variable;
	 vector absty)
{
    INIT_TRACE;

    vector env = resolvemissedray(dir, Time, mask);

    if (maxdist == .0)
	eval = env;
    else if (style)
	{
	    if (style == 1)
		if (doabs)
		    {
			if(TRACE_ABSRP)
			    eval = hitCf * exp(-raylength * absty);
		    }
		else
		    {
			if (TRACE_GATHER) eval = hitCf;
			else eval = env;
		    }
	    else
		if (TRACE_OCCL) eval = .0;
		else eval = env;

	}

    return max(.0, eval);
}

// Cone case
vector
raytrace(vector p, dir;
	 float angle, maxdist;
	 int samples, oblend, style;
	 string raystyle, scope, variable;
	 vector absty)
{
    INIT_TRACE;

    vector env = resolvemissedray(dir, Time, mask, "angle", angle);

    if (maxdist == .0)
	eval = env;
    else if (style)
	if (style == 1)
	    {
		vector tmp = .0;
		vector raydir;

		if (doabs)
		    {
			gather(p, dir,
			       TRACE_FLAGS(samples),
			       "angle", angle,
			       variable, hitCf,
			       "variancevar", variable,
			       "ray:length", raylength)
			    {
				tmp += hitCf * exp(-raylength * absty);
			    }
		    }
		else
		    {
			gather(p, dir,
			       TRACE_FLAGS(samples),
			       "angle", angle,
			       variable, hitCf,
			       "variancevar", variable,
			       "ray:direction", raydir)
			    {
				tmp += hitCf;
			    }
			else
			    {
				tmp += resolvemissedray(raydir, Time, mask);
			    }
		    }

		eval = tmp / samples;
	    }
	else
	    eval = 2.0 * occlusion(p, dir,
				   TRACE_FLAGS(samples),
				   "background", 1.)
		* env;

    return max(.0, eval);
}


// Use BSDF
#define SAMPLE_BSDF sample_bsdf(f, v, dir, brdf, pdf, type, sx, sy, mask)
#define CONTRIBUTE eval += pdf * brdf * tmp; summ += pdf
#define AVERAGE eval /= summ

#define VARIANCEAA if (dorayvariance)					\
	{ float lum = luminance(eval) / summ;				\
	    if (isgamma) lum = sqrt(lum);				\
	    int _i1 = _i + 1;						\
	    if (_i1 >= minraysamples)					\
		{ int samplesize;					\
		    float mean;						\
		    float newvar = variance(lum - prevlum,		\
					    mean, samplesize);		\
		    var = (var * _i + newvar) / (_i1);			\
		    if (var <= variance*variance) break; }		\
	    prevlum = lum; }

#define FINALIZE_SAMPLING CONTRIBUTE; VARIANCEAA; END_LOOP; AVERAGE

vector
raytrace(bsdf f;
	 vector p, v;
	 float maxdist;
	 int sid, samples, oblend, style;
	 string raystyle, scope, variable;
	 int dorayvariance, minraysamples;
	 float variance;
	 int isgamma;
	 vector absty)
{
    INIT_TRACE;

    int type;
    vector dir, brdf, tmp;
    float pdf;
    float summ = .0;
    float prevlum = .0;
    float var = .0;

    if (maxdist == .0)
	{
	    START_SAMPLING("nextpixel");
	    SAMPLE_BSDF;

	    tmp = resolvemissedray(dir, Time, mask);

	    FINALIZE_SAMPLING;
	}
    else if (style)
	if (style == 1)
	    {
		START_SAMPLING("nextpixel");
		SAMPLE_BSDF;

		if (doabs)
		    {
			if(TRACE_ABSRP)
			    tmp = hitCf * exp(-raylength * absty);
		    }
		else
		    {
			if (TRACE_GATHER)
			    tmp = hitCf;
			else
			    tmp = resolvemissedray(dir, Time, mask);
		    }

		FINALIZE_SAMPLING;
	    }
	else
	    {
		START_SAMPLING("nextpixel");
		SAMPLE_BSDF;

		if (TRACE_OCCL)
		    tmp = .0;
		else
		    tmp = resolvemissedray(dir, Time, mask);

		FINALIZE_SAMPLING;
	    }

    return eval;
}


#define INIT_ABSRP float raylength = 1.;	\
    string sfilter = "closest";			\
    string raystyle = "refract";		\
    vector eval = .0

// volume absorption
// Delta case
vector
absorption(vector p, dir, kabs;
	   float maxdist;
	   string scope)
{
    INIT_ABSRP;

    if(trace(p, dir, Time,
	     TRACE_FLAGS(1),
	     "ray:length", raylength))
	eval = exp(-raylength * kabs);

    return eval;
}

// Cone case
vector
absorption(vector p, dir, kabs;
	   float maxdist, angle;
	   int samples;
	   string scope)
{
    INIT_ABSRP;

    gather(p, dir,
	   TRACE_FLAGS(samples),
	   "angle", angle,
	   "ray:length", raylength)
	{
	    eval += exp(-raylength * kabs);
	}

    return eval / samples;
}


// Shadow absorption
vector
shadowabs(vector p, absty)
{
    float bias;
    renderstate("renderer:raybias", bias);

    float raylength = rayhittest(p, -I, bias);

    if (raylength < .0)
        raylength = length(I);

    return exp(-raylength * absty);
}


// Dielectric fresnel at normal incidence
float
fresnel_ni(float eta)
{
    float eval = (eta - 1.)/(eta + 1.);
    return eval * eval;
}


// Dielectric fresnel (Schlick approximation) albedo
void
fresnel_albedo(float eta;
	       export float ar, at)
{
    ar = fresnel_ni(eta) * 0.8333333333 + 0.1666666666;
    at = 1. - ar;
}


// Simple energy distribution
void
edist(float alb, eta;
      int enableDFS, enableSPC, enableTRN, enableSSS;
      float wdfs, wspc, wtrn, wsss;
      export float kdfs, kspc, ktrn, ksss)
{
    float summ,
	wgloss = 0,
	_spc = enableSPC * wspc,
	_trn = enableTRN * wtrn,
	_dfs = enableDFS * wdfs,
	_sss = enableSSS * wsss;

    // Reflection and refraction are depened on incident ray
    wgloss = max(enableSPC, enableTRN) * max(wspc, wtrn);

    _spc /= wgloss;
    _trn /= wgloss;

    // Conserve fresnel attenuated energy
    if(enableSPC || enableTRN)
	{
	    float ar, at;
	    fresnel_albedo(eta, ar, at);
	    _dfs /= at;
	}

    // Diffuse and SSS are uniform on hemisphere
    summ = _dfs + _sss;

    if (enableTRN)
	{
	    kspc = alb * _spc;
	    ktrn = alb * _trn;
	}
    else
	{
	    summ = max(summ, wgloss);

	    _spc /= summ;
	    kspc = alb * _spc;
	}

    _dfs /= summ;
    _sss /= summ;
    kdfs = alb * _dfs;
    ksss = alb * _sss;
}


// Calculates the point where refracted ray goes out
// and the amount of absorption
void
thinP(vector p, i, nbN, nfN;
      float absty, eta, thickness;
      export vector newP;
      export float absrp)
{
    vector tdir = refract(i, nfN, eta);
    float len = thickness / abs(dot(tdir, nbN));
    newP = p + tdir * abs(len);
    absrp = exp(-absty * len);
}


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

    START_SAMPLING("nextpixel");
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

    START_SAMPLING("nextpixel");

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
	    vector hitCf = illum_surface(hitP, hitN, eta, sid,
					depth, depthimp);
	    eval += weight * hitCf * bssrdf->eval(Xr, n);
	    pdf += weight;
	}

    END_LOOP;

    return eval / pdf;
}


// Invert hue of given RGB
vector
invert_hue(vector color)
{
    vector hsv = rgbtohsv(color);
    hsv.x = (hsv.x + 0.5) % 1;
    return hsvtorgb(hsv);
}


// The main surface routine
void
physurface(int conductor;
	   int thin;
	   float thickness;	// Thickness of thin wall if thin
	   float alb;		// Albedo
	   vector clrsurf;	// Surface color
	   vector ior;		// Tuple of refraction indicies
	   int enableDFS, enableSPC, enableTRN, enableSSS;
	   float weightDFS, weightSPC, weightTRN, weightSSS;
	   float roughDFS, roughSPC;
	   float anisobias;
	   vector clrSSS;	// Scattering coefficient
	   vector absty;	// Absorption coefficient
	   float g;		// Scattring phase
	   int dispersion;
	   int styleSPC, styleTRN;   // How to perform reflection/refraction
	   int oblendSPC, oblendTRN; // Opacity blending
	   int squality;	// Sampling quality (see sampling_quality)
	   int tsamples;	// Number of ray-tracing samples
	   int vsamples;	// Number of single scattering samples
	   int ssamples;	// Number of multiple scattering samples
	   int empty;
	   int useF;		// Use BSDF to compute reflection/refraction
	   int accurateabs;	// Accurate glossy absorption
	   float depthimp;	// Depth importance
	   float maxdist;	// Maximum tracing distance
	   vector p;		// Point position
	   vector n;		// Point normal
	   vector ii;		// Incident ray from eye to point
	   vector tangent;	// Surface derivative
	   vector sellmeierB, sellmeierC; // Sellmeier's coefficients
	   string gvarSPC, gvarTRN;	  // Gather variables
	   string sscope;		  // Multiple scattering object scope
	   export vector beauty;
	   export vector opacity;
	   export bsdf f;
	   export vector all[])
{
    beauty = .0;
    opacity = 1.;
    vector fullDFS = .0;
    vector fullSPC = .0;
    vector fullTRN = .0;
    vector fullSSS = .0;

    
    string renderengine;
    renderstate("renderer:renderengine", renderengine);
    int sid = renderengine == "micropoly" ? newsampler() : SID;


    float
	kDFS = .0,
	kSPC = .0,
	kTRN = .0,
	kSSS = .0;

    // Fresnel reflectance
    float
	fr = 1.,
	ft = 1.;

    // IOR
    float etat = max(1.0000016, ior.x);

    // Dispersion
    vector dtint = 1.;
    if (dispersion)
	{
	    float wl = samplewl(rand(sid));
	    etat = sellmeier(wl, sellmeierB, sellmeierC);
	    dtint = wl2rgb(wl);
	}

    float
	etai = ior.y,
	etak = ior.z,
	eta = etai/etat;

    // Geometry atenuation factor
    float
	gafmask = 1.,
	gafrefr = 1.;

    // Specular BRDF parameter
    // With raytrace | micropoly engine
    // light specular still evaluates
    float sigma = max(roughSPC, SMOOTH_THRESHOLD);

    // Treat surface as perfect smooth
    // if specular rougness is below this value
    int smooth = roughSPC <= SMOOTH_THRESHOLD;

    // When reflections/refraction are anisotropic
    int isanisotropic = (anisobias != .0) && !smooth;

    // Anisotropy roughness
    float
	sigmaU = sigma,
	sigmaV = sigma;
    
    // Ray-tracing spread angle
    float angle = .0;

    // Absorption coefficient, scalar
    float kabs = max(absty);

    vector ni = normalize(ii);
    //Vector from point to viewer
    vector v = -ni;

    // Refraction outgoing ray point for thin plate
    vector pTRN = p;

    vector nfN = frontface(n, ni);
    vector nbN = -nfN;

    vector rdir = reflect(ni, nfN);
    vector tdir = thin ? ni : refract(ni, n, eta);

    vector
	tU = tangent,
	tV = cross(tU, n);

    // Normalized surface color
    vector _clrsurf = clrsurf / ALONE_VEC(clrsurf);
    vector clrDFS = _clrsurf;
    vector clrSPC = 1.0;

    // Normalized refraction color
    vector clrTRN = absty / ALONE(kabs);

    // Exponent inverts color. Protect from this.
    vector _absty = invert_hue(absty);

    // Amount of refraction absorption if thin
    float _absTRN = 1.;
    // and else
    vector absTRN = 1.;

    // Amount of reflection absorption
    // for internal reflection
    vector absSPC = 1.;

    vector factorDFS = .0;
    vector factorSPC = .0;
    vector factorSSS = .0;
    vector factorTRN = .0;

    // Ray-traced values
    vector traceTRN = .0;
    vector traceSPC = .0;

    // Shadow amount
    vector sh = .0;

    // Absorption direction
    vector absdir = tdir;

    // Single scattering eval
    vector singlescattering = .0;

    // Ray is directed inside object
    int	enter = dot(n, ni) < .0;

    int depth = getraylevel() + getglobalraylevel();

    // tracing samples
    float variance;
    renderstate("object:variance", variance);
    int dorayvariance;
    renderstate("object:dorayvariance", dorayvariance);
    string colorspace;
    renderstate("renderer:colorspace", colorspace);
    int isgamma = colorspace == "gamma";
    int maxraysamples, minraysamples;
    renderstate("light:maxraysamples", maxraysamples);
    renderstate("light:minraysamples", minraysamples);

    // Copied variables starts with "_"
    int _tsamples = tsamples;
    if (!smooth)
	if (useF && !squality && dorayvariance)
	    // _tsamples as max ray samples
	    _tsamples = maxraysamples;
	else
	    {
		if (squality == 0)
		    _tsamples = floor(lerp(minraysamples, maxraysamples,
					   max(sigma / MAX_ROUGH,
					       MAX_ROUGH)));
		else if (squality == 1) _tsamples = minraysamples;
		else if (squality == 2) _tsamples = maxraysamples;
	    }

    int _vsamples = vsamples;
    int _ssamples = ssamples;

    // Is the total internal reflection case
    int internal = rdir == tdir;

    int solid = !thin;
    int thick = thin && thickness > .0;

    // Recompute number of samples by the depth importance
    if (depth && depthimp < 1.)
	{
	    float factor = pow(depthimp, depth);

	    _tsamples = FLOOR_ALONE(tsamples * factor);
	    _ssamples = FLOOR_ALONE(ssamples * factor);
	    _vsamples = FLOOR_ALONE(vsamples * factor);
	}

    // Ray-tracing scope
    string scopeSPC = "scope:default";
    string scopeTRN = "scope:default";

    bsdf f_DFS = diffuse(nfN, roughDFS);
    bsdf f_SPC = specular(rdir);
    bsdf f_TRN = specular(tdir);
    bsdf f_SSS = diffuse(nbN);
    bsdf f_VOL = g == .0 ? isotropic() : henyeygreenstein(g);


    // When specular is allowed
    // In case of total internal reflection
    // specular is computed in refraction routine
    int allowSPC =
	enableSPC	&&
	!internal;

    // When allow transmission
    int allowTRN =
	!conductor	&&
	!enableDFS	&&
	enableTRN;

    // When SSS is allowed
    int allowSSS =
	enableSSS	&&
	solid;

    int allowmultisss =
	allowSSS	&&
	_ssamples;

    int allowsinglesss =
	allowSSS	&&
	_vsamples;

    int translucent =
	enableSSS	&&
	thin;

    // Single scattering
    RayMarcher singlesss;
    singlesss->init(_absty, f_VOL, sid, _vsamples, depth, depthimp);

    // disable separate absorption and single scattering
    // for Raytrace/Micropoly renderers
    int doAbsTRN = enter || internal;

    int oAbs =
	kabs > .0	&&
	(renderengine == "micropoly" ||
	 renderengine == "raytrace");

    int oAbsDeltaTRN = oAbs && doAbsTRN && styleTRN;

    // internal reflection
    int oAbsDeltaSPC = oAbs && !doAbsTRN && styleSPC;

    // glossy absorption
    int oAbsGloss = !smooth && accurateabs && oAbs;

    int oAbsGlossTRN = oAbsGloss && oAbsDeltaTRN;
    int oAbsGlossSPC = oAbsGloss && oAbsDeltaSPC;

    // Get tracing masks for non-gather tracing
    if(!renderstate("object:reflectmask", scopeSPC))
	scopeSPC = "scope:default";
    if(!renderstate("object:refracttmask", scopeTRN))
	scopeTRN = "scope:default";

    if (empty && solid)
	{
	    if (enter || internal)
		scopeTRN = "scope:self";
	    if (!enter)
		scopeSPC = "scope:self";
	}

    // The main coefficients
    if(conductor)
	{
	    kSPC = alb;
	    fr = cfresnel(v, nfN, eta, etak);
	    clrSPC = _clrsurf;
	}
    else
	{
	    edist(alb, eta,
		  enableDFS, enableSPC, enableTRN, enableSSS,
		  weightDFS, weightSPC, weightTRN, weightSSS,
		  kDFS, kSPC, kTRN, kSSS);

	    // This also affects the shadow
	    // so use frontfaced normal
	    fresnel(ni, nfN, eta, fr, ft);
	}

    // Diffuse
    if (enableSPC)
	f_DFS = cvex_bsdf("fresneldiffuse_eval",
			  "fresneldiffuse_sample",
			  "label", "diffuse",
			  "N", nfN,
			  "eta", eta);

    // BSDFs of reflections and refractions
    if(smooth)
	{
	    f_SPC = specular(rdir);
	    f_TRN = specular(tdir);
	}
    else
	{
	    angle = atan(sigma);

	    if(isanisotropic)
		{
		    anisorough(sigma, anisobias, sigmaU, sigmaV);

		    float beta = anisorough(v, tU, tV, sigmaU, sigmaV);

		    gafmask = gaf(dot(v, nfN), beta);
		    gafrefr = gaf(dot(tdir, nbN), beta);

		    f_SPC = cvex_bsdf("phy_aniso_eval",
				      "phy_aniso_sample",
				      "label", "reflect",
				      "n", nfN,
				      "sigmau", sigmaU,
				      "sigmav", sigmaV,
				      "tu", tU,
				      "tv", tV);

		    f_TRN = cvex_bsdf("phy_aniso_trans_eval",
				      "phy_aniso_trans_sample",
				      "label", "refract",
				      "n", n,
				      "sigmau", sigmaU,
				      "sigmav", sigmaV,
				      "tu", tU,
				      "tv", tV,
				      "eta", thin ? 1. : eta);
		}
	    else
		{
		    gafmask = gaf(dot(v, nfN), sigma);
		    gafrefr = gaf(dot(tdir, nbN), sigma);

		    f_SPC = cvex_bsdf("phy_eval",
		    		      "phy_sample",
		    		      "label", "reflect",
		    		      "n", nfN,
		    		      "sigma", sigma,
				      "tu", tU,
				      "tv", tV);

		    f_TRN = cvex_bsdf("phy_trans_eval",
				      "phy_trans_sample",
				      "label", "refract",
				      "n", n,
				      "sigma", sigma,
				      "tu", tU,
				      "tv", tV,
				      "eta", thin ? 1. : eta);
		}
	}


    // Ray-tracing specular
    if (allowSPC && styleSPC)
	if (smooth)
	    traceSPC = raytrace(p, rdir,
				maxdist,
				oblendSPC, styleSPC,
				"reflect", scopeSPC, gvarSPC,
				oAbsDeltaSPC ? _absty : .0);
	else if (useF)
	    traceSPC = raytrace(f_SPC,
				p, v,
				maxdist,
				sid, _tsamples, oblendSPC, styleSPC,
				"reflect", scopeSPC, gvarSPC,
				dorayvariance, minraysamples,
				variance, isgamma,
				oAbsGlossSPC ? _absty : .0);
	else
	    traceSPC = raytrace(p, rdir,
				angle, maxdist,
				_tsamples, oblendSPC, styleSPC,
				"reflect", scopeSPC, gvarSPC,
				oAbsGlossSPC ? _absty : .0);

    // Refraction
    if (allowTRN && styleTRN)
	{
	    int do_trace = 1;

	    if (thin)
		if (thick)
		    thinP(p, ni, nbN, nfN,
			  kabs, eta, thickness,
			  pTRN, _absTRN);
		else
		    {
			if (smooth)
			    {
				opacity = 1. - kTRN * ft * clrTRN;

				// Don't pathtrace refraction
				// because already used opacity
				f_TRN *= .0;
				do_trace = 0;
			    }
		    }
	    else if (kabs > .0)
		{
		    vector abstmp = 1.;
		    vector tmpsss = .0;

		    if (!oAbs)
			{
			    if (accurateabs && !smooth)
				abstmp = absorption(p, absdir,
						    _absty,
						    maxdist, angle,
						    _tsamples,
						    scopeTRN);
			    else
				abstmp = absorption(p, absdir,
						    _absty,
						    maxdist,
						    scopeTRN);
			}

		    if (allowsinglesss)
			{
			    float raylength = .0;

			    if (trace(p, absdir, Time,
				      "samplefilter", "closest",
				      "raystyle", "refract",
				      "scope", scopeTRN,
				      "maxdist", maxdist,
				      "ray:length", raylength))
				{
				    vector clr = clrSSS / ALONE_VEC(clrSSS);
				    tmpsss = singlesss->eval(p, absdir,
							     raylength)
					* clr;
				}
			}

		    if (enter || internal)
			{
			    absTRN = abstmp;
			    singlescattering = tmpsss;
			}
		    else if(enableSPC)
			{
			    absSPC = abstmp;
			    singlescattering = tmpsss;
			}
		}

	    if (do_trace)
		if (smooth)
		    traceTRN = raytrace(pTRN, tdir,
					maxdist,
					oblendTRN, styleTRN,
					"refract", scopeTRN, gvarTRN,
					oAbsDeltaTRN ? _absty : .0);
		else if (useF)
		    traceTRN = raytrace(f_TRN,
					pTRN, v,
					maxdist,
					sid, _tsamples, oblendTRN, styleTRN,
					"refract", scopeTRN, gvarTRN,
					dorayvariance, minraysamples,
					variance, isgamma,
					oAbsGlossTRN ? _absty : .0);
		else
		    traceTRN = raytrace(pTRN, tdir,
					angle, maxdist,
					_tsamples, oblendTRN, styleTRN,
					"refract", scopeTRN, gvarTRN,
					oAbsGlossTRN ? _absty : 0);
	}


    // The base factors
    factorDFS = clrDFS * kDFS;
    factorSPC = clrSPC * kSPC * fr * gafmask;
    factorTRN = kTRN * ft * gafrefr * dtint;
    factorSSS = kSSS;

    // Translucency
    if (thin)
	{
	    factorSSS *= clrSSS / ALONE_VEC(clrSSS);
	    factorTRN *= clrTRN * _absTRN;
	}

    // Lighting
    illum_surface(p, pTRN, p + nbN * thickness,
    		  nbN, nfN,
    		  v,
    		  thick,
    		  thickness,
    		  depthimp,
    		  sid,
    		  depth,
    		  enableDFS, allowSPC, allowTRN, translucent,
    		  f_DFS, f_SPC, f_TRN, f_SSS,
    		  factorDFS, factorSPC, factorTRN, factorSSS,
    		  fullDFS, fullSPC, fullTRN, fullSSS);

    // Compute multiple scattering
    if (allowmultisss)
	fullSSS = raySSS(p, n,
			 eta, g,
			 _absty, clrSSS,
			 _ssamples, sid,
			 sscope,
			 depth, depthimp)
	    * factorSSS;

    if (enableTRN)
	{
	    if (allowsinglesss)
		fullSSS += singlescattering * factorSSS;

	    fullDFS *= factorTRN;
	    fullSSS *= factorTRN;

	    f_DFS *= factorTRN;
	    f_SSS *= factorTRN;
	}

    
    if (allowTRN)
	{
	    fullTRN += traceTRN * absTRN * factorTRN;
	    f_TRN *= factorTRN * absTRN;
	}
    else
	{
	    fullTRN = .0;
	    f_TRN *= .0;
	}

    
    fullSPC += traceSPC * absSPC * factorSPC;

    // PBR BSDF's
    f_DFS *= factorDFS;
    f_SPC *= factorSPC * absSPC;
    f_SSS *= factorSSS;

    // Transparent shadow
    if (isshadowray() && allowTRN)
	{
	    sh = ft * kTRN;

	    if (thin)
		{
		    if(thick)
			sh *= exp(-kabs * thickness / dot(ni, nbN));

		    sh *= clrTRN;
		}
	    else if (!enter)
		sh *= shadowabs(p, _absty);

	    opacity = 1. - sh;
	}

    beauty = fullDFS + fullSPC + fullTRN + fullSSS;
    all = array(fullDFS, fullSPC, fullTRN, fullSSS);
    f = f_DFS + f_SPC + f_TRN + f_SSS;
}


#endif // __physurface__
