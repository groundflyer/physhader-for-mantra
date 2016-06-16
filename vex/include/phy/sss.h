// This may look like -*- C -*- code, but it is really Houdini Vex
//
//	sss.h - subsurface scattering implementation,
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


#ifndef __physss__
#define __physss__


#include <expsampler.h>
#include <phy/utils.h>


// Irradiance for sss
vector
illum_surface(vector p, n;
	      float eta;
	      int sid;
	      int depth;
	      float depthimp;
	      int doshadow;
	      string lightmask)
{
    vector eval = .0;

    foreach (int lid; getlights("lightmask", lightmask))
	{
	    int mask, samples = 1;
	    if (setcurrentlight(lid)) {
		int isarealight = 0;
		renderstate("light:arealight", isarealight);
		if (isarealight) {
		    renderstate("light:maxraysamples", samples);
		    if (depth && depthimp != 1.)
			samples = FLOOR_ALONE(samples * pow(depthimp,depth)); } }
    
	    vector accum = .0;

	    START_SAMPLING("nextpixel");
	    SET_SAMPLE;

	    float scale;
	    vector lp, leval;
	    mask = sample_light(lid, p, sample,
				Time, lp, leval, scale);

	    vector cl = leval * scale / ALONE_VEC(leval);

	    vector l = lp - p;

	    if (doshadow)
		cl *= shadow_light(lid, p, l, Time,
				   "noantialiasing", 1,
				   "N", n,
				   "SID", sid);

	    if(mask & bouncemask("sss"))
		{
		    float kr, kt;
		    l = normalize(-l);
		    fresnel(l, n, eta, kr, kt);
		    accum += cl * kt * max(dot(n, -l), .0);
		}

	    END_LOOP; 			// SAMPLING

	    eval += accum/samples;

	}

    return eval;
}


// builds a transform matrix randomly rotated along axiz z
matrix3
rbasis(float sx; vector z)
{
    vector x, y;
    vector u = normalize(vector(rand(sx))-0.5);
    makebasis(x, y, z, u);
    return set(x, y, z);
}


// evaluates BSSRDF approximate reflectance profile with given albedo and radius
vector
reflectance_profile(vector alb;
		    float r)
{
    vector tmp = -r/alb;
    return (exp(tmp) + exp(tmp/3))/(1.+r);
}


// returns relative position for sampling
vector
sss_sample_pos(float sx, sy, alb, radius)
{
    float theta = 2. * PI * sx;
    float tmp = -sy*radius/alb;
    float z = 0.25*exp(tmp)+0.75*exp(tmp/3);
    return set(sy*cos(theta), sy*sin(theta), z);
}


// computes subsurface scattering
vector
sss_multi(vector p;
	  vector n;
	  vector alb;
	  float eta;
	  int sid;
	  string scope;
	  int doshadow;
	  float curvature;
	  string lightmask;
	  int depth;
	  float depthimp;
	  VarianceSampler vsampler)
{
    float falb = max(alb);
    float radius = 7.50184474 * pow(falb, 0.78677001);

    expsampler samp;
    samp->init(alb, radius);

    vector eval = .0;
    float pdf = .0;

    int samples = vsampler.maxraysamples;

    START_SAMPLING("decorrelate");

    vector randn = normalize(vector(rand(sx))-0.5);
    vector sn = lerp(n, randn, curvature);
    vector pt = radius * sss_sample_pos(sx, sy, falb, radius);
    vector dir = set(.0, .0, -1.);
    matrix3 basis = rbasis(sy, sn);
    pt = p + ptransform(pt, basis);
    dir = ntransform(dir, basis);

    vector hitP, hitN;

    if (trace(pt, dir, Time,
              "SID", sid,
              "bias", 0,
              "pipeline", "displacement",
              "scope", scope,
              "samplefilter", "closest",
              "maxdist", radius,
              "P", hitP,
              "N", hitN))
        {
            float r = distance(p, hitP);
            vector irr = illum_surface(hitP, normalize(hitN),
				       eta,
				       sid, depth, depthimp,
				       doshadow, lightmask);

            float sval = sy * samp.k3 * samp.max_rand;
            int icomp = sampleExpComp(sval, samp.k1, samp.k2);

            vector evalR = reflectance_profile(alb, r);
            float weight = 1. / exp(-samp.ext_sort[icomp] * r) * samp.ext_sort[icomp];

            eval += weight * evalR * irr;
            pdf += weight;
        }

    if (vsampler->stop_by_variance(max(eval)/pdf, _i))
	    break;

    END_LOOP;

    return eval / pdf;
}


// volume direct lighting
vector
illum_volume(vector p, v;
	     bsdf f;
	     int sid;
	     int depth;
	     float depthimp;
	     int shadow;
	     string lightmask)
{
    vector eval = .0;

    vector l, cl;
    foreach (int lid; getlights("lightmask", lightmask))
	{
	    int mask, samples = 1;
	    if (setcurrentlight(lid)) {
		int isarealight = 0;
		renderstate("light:arealight", isarealight);
		if (isarealight) {
		    renderstate("light:maxraysamples", samples);
		    if (depth && depthimp != 1.)
			samples = FLOOR_ALONE(samples * pow(depthimp,depth)); } }

	    vector accum = .0;
	    float pdf = 0;

	    START_SAMPLING("nextpixel");
	    SET_SAMPLE;

	    SAMPLE_LIGHT(p, v);

	    float weight = 0;
	    cl *= eval_bsdf(f, v, l, weight, PBR_VOLUME_MASK);

	    pdf += weight;
	    accum += cl * weight;

	    END_LOOP; 	// SAMPLING

	    eval += accum / pdf;
	}

    return eval;
}


// Constant density stachaostic raymarcher
struct RayMarcher
{
    vector ca; 			// absorption coefficient
    bsdf f;			// phase function
    int sid;
    int samples;		// number of samples
    int depth;
    float depthimp;
    int doshadow;
    string lightmask;

    VarianceSampler vsampler;


    vector
    eval(vector p, v;
	 float raylength)
    {
	vector accum = .0;
	int counter = 1;

	expsampler samp;
	samp->init(ca, raylength);

	int samples = vsampler.maxraysamples;

	START_SAMPLING("decorrelate");

	vector cl;
	float spo = samp->sample(cl, sx);

	vector pp = p + v * spo;

	cl *= illum_volume(pp, v, f, sid, depth, depthimp, doshadow, lightmask);

	accum += cl;

	if (vsampler->stop_by_variance(max(accum)/(counter++), _i))
	    break;
	
	END_LOOP;

	return accum / counter;
    }
};


#endif	// _physss__
