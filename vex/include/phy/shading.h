// This may look like -*- C -*- code, but it is really Houdini Vex
//
//	shading.h - Surface shading routines for raytrace/micropoly engine.
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


#ifndef __phy_shading__
#define __phy_shading__


#include <math.h>
#include <pbr.h>
#include <phy/utils.h>


// Treat surface smooth if roughness below this value
#define SMOOTH_THRESHOLD	.00005


// Illuminance loop
#define START_ILLUMINANCE				  \
    vector l, cl;					  \
    foreach (int lid; getlights()) {			  \
    int mask, samples = 1;						\
    if (setcurrentlight(lid)) {						\
    int isarealight = 0;						\
    renderstate("light:arealight", isarealight);			\
    if (isarealight) {							\
    renderstate("light:maxraysamples", samples);			\
    if (depth && depthimp != 1.)					\
	samples = FLOOR_ALONE(samples * pow(depthimp,depth)); } }


// Eval bsdf with current sample
#define EVAL_BSDF(x)	eval_bsdf(x, v, l, mask)


// Sample light with given position and normal
#define SAMPLE_LIGHT(PP, NN)			\
    { mask = sample_light(lid, sid,		\
			  PP, NN, sample,	\
			  cl, l);		\
	l = normalize(l); }


// Wrapper for sample_light and shadow_light
int
sample_light(int lid, sid;
	     vector p, n, sample;
	     export vector cl, l)
{
    vector lp, eval;
    float scale;

    int mask = sample_light(lid, p, sample,
			    Time, lp, eval, scale);

    cl = eval * scale / ALONE_VEC(eval);
    l = lp - p;
    cl *= shadow_light(lid, p, l, Time,
		       "N", n,
		       "SID", sid);

    return mask;
}


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

    START_SAMPLING;
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


// Conductor fresnel
float
cfresnel(vector v, n;
	 float eta, k)
{
    float
	cosi = dot(v, n),
	tmp = eta * eta + k * k,
	cosi2 = cosi * cosi,
	tmp1 = tmp * cosi2,
	tmp2 = 2 * eta * cosi,
	par = (tmp1 - tmp2 + 1.0) / (tmp1 + tmp2 + 1.0),
	per = (tmp - tmp2 + cosi2) / (tmp + tmp2 + cosi2);
    return (par + per) / 2.0;
}


// Volume absorption
// Delta case
vector
absorption(vector p, dir, kabs;
	   float maxdist;
	   string scope)
{
    vector eval = .0;
    float raylength = .0;
    if(trace(p, dir, Time,
	     "samplefilter", "closest",
	     "scope", scope,
	     "samples", 1,
	     "maxdist", maxdist,
	     "raystyle", "refract",
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
    vector eval = .0;
    float raylength = .0;
    gather(p, dir,
	   "samplefilter", "closest",
	   "scope", scope,
	   "samples", samples,
	   "maxdist", maxdist,
	   "angle", angle,
	   "raystyle", "refract",
	   "ray:length", raylength)
	{
	    eval += exp(-raylength * kabs);
	}

    return eval / samples;
}


// Trace reflection/refraction
// Style choices:
//	1 - Gather
//	2 - Occlusion
// Delta function case
vector
raytrace(vector p, dir;
	 float maxdist;
	 int oblend, style;
	 string raystyle, scope, gathervar)
{
    vector hitCf = .0;

    int mask = PBR_REFRACT_MASK;

    if (raystyle == "reflect")
	mask = PBR_REFLECT_MASK;

    vector eval = resolvemissedray(dir, Time, mask);

    string sfilter = oblend ? "opacity" : "closest";

    if (maxdist != .0)
	if (style == 1)
	    {
		if (trace(p, dir, Time,
			  "scope", scope,
			  "samplefilter", sfilter,
			  "maxdist", maxdist,
			  "samples", 1,
			  "raystyle", raystyle,
			  gathervar, hitCf))
		    eval = hitCf;
	    }
	    else
		if (trace(p, dir, Time,
			  "scope", scope,
			  "samplefilter", sfilter,
			  "maxdist", maxdist,
			  "samples", 1,
			  "raystyle", raystyle))
		    eval = .0;

    return eval;
}

// Cone case
vector
raytrace(vector p, dir;
	 float angle, maxdist;
	 int samples, oblend, style;
	 string raystyle, scope, gathervar)
{
    vector hitCf = .0;
    vector eval = .0;
    int mask = PBR_REFRACT_MASK;

    if (raystyle == "reflect")
	mask = PBR_REFLECT_MASK;

    vector env = resolvemissedray(dir, Time, mask, "angle", angle);

    string sfilter = oblend ? "opacity" : "closest";

    if (maxdist == .0)
	eval = env;
    else if (style)
	if (style == 1)
	    {
		vector tmp = .0;
		vector raydir;

		gather(p, dir,
		       "scope", scope,
		       "samplefilter", sfilter,
		       "angle", angle,
		       "maxdist", maxdist,
		       "samples", samples,
		       "raystyle", raystyle,
		       gathervar, hitCf,
		       "variancevar", gathervar,
		       "ray:direction", raydir)
		    {
			tmp += hitCf;
		    }
		else
		    {
			tmp += resolvemissedray(raydir, Time, mask);
		    }

		eval = tmp / samples;
	    }
	else
	    eval = 2.0 * occlusion(p, dir,
				   "samplefilter", sfilter,
				   "angle", angle,
				   "maxdist", maxdist,
				   "samples", samples,
				   "scope", scope,
				   "background", 1.)
		* env;

    return max(.0, eval);
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
    float ar, at, summ,
	wgloss = 0,
	_spc = enableSPC * wspc,
	_trn = enableTRN * wtrn,
	_dfs = enableDFS * wdfs,
	_sss = enableSSS * wsss;

    wgloss = max(enableSPC, enableTRN) * max(wspc, wtrn);

    _spc /= wgloss;
    _trn /= wgloss;

    // Conserve fresnel attenuated energy
    if(enableSPC || enableTRN)
	{
	    fresnel_albedo(eta, ar, at);
	    _dfs /= at;
	}

    summ = _dfs + _sss;

    if (enableTRN)
	{
	    kspc = alb * _spc;
	    ktrn = alb * _trn;
	}
    else
	{
	    summ += _spc;
	    kspc = alb * _spc / summ;
	}

    kdfs = alb * _dfs / summ;
    ksss = alb * _sss / summ;
}


// Calculates the point where refracted ray goes out
// and the amount of absorption
void
thinP(vector p, i, nfN;
      float absty, eta, thickness;
      export vector newP;
      export float absrp)
{
    vector tdir = refract(i, nfN, eta);
    float len = thickness / dot(tdir, -nfN);
    newP = p + tdir * len;
    absrp = exp(-absty * len);
}


#endif // __phy_shading__
