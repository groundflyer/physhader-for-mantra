// This may look like -*- C -*- code, but it is really Houdini Vex
//
//	utils.h - Common utilities for phyShader.
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


#ifndef __phy_utils__
#define __phy_utils__


// At Least One
#define ALONE(x)	max(x, 1)
#define FLOOR_ALONE(x)	ALONE(floor(x))
#define ALONE_VEC(x)	ALONE(max(x))

#define START_SAMPLING(MODE)					\
    float sx, sy;						\
    for (int _i = 0; _i < samples; ++_i) {			\
    nextsample(sid, sx, sy, "mode", MODE)


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


// Sample light with given position and normal
#define SAMPLE_LIGHT(PP, NN)			\
    mask = sample_light(lid, sid, shadow,	\
			PP, NN, sample,		\
			cl, l);			\
    l = normalize(l);

#define SET_SAMPLE	vector sample = set(sx, sy, .0);

#define END_LOOP	}

#define FIT01(X) fit(X, min(X), max(X), .0, 1.)


// Wrapper for sample_light and shadow_light
int
sample_light(int lid, sid;
	     int shadow;
	     vector p, n, sample;
	     export vector cl, l)
{
    vector lp, eval;
    float scale;

    int mask = sample_light(lid, p, sample,
			    Time, lp, eval, scale, "N", n);

    cl = eval * scale / ALONE_VEC(eval);
    l = lp - p;

    if (shadow)
	cl *= shadow_light(lid, p, l, Time,
			   "N", n,
			   "SID", sid);

    return mask;
}


// Test for point is in object
int
inobject(vector p;
         string scope;
	 vector test_dir)
{
    vector hitN;
    int hit = trace(p, test_dir, Time, "scope", scope, "N", hitN);
    return hit && (dot(test_dir, hitN) > .0);
}


// Invert hue of given RGB
vector
invert_hue(vector color)
{
    vector hsv = rgbtohsv(color);
    hsv.x = (hsv.x + 0.5) % 1;
    return hsvtorgb(hsv);
}


// Variance anti-aliasing helper function
int
stop_by_variance(float _lum, prevlum, variance;
		 int isgamma, isample, minraysamples;
		 export float var)
{
    int i1 = isample + 1;

    if (i1 >= minraysamples)
	{
	    float lum = _lum;
	    if (isgamma) lum = sqrt(lum);

	    int samplesize;
	    float mean;
	    float newvar = variance(lum - prevlum,
				    mean, samplesize);
	    var = (var * isample + newvar) / i1;

	    if (var <= variance*variance)
		return 1;
	}
    return 0;
}


struct VarianceSampler
{
    int dorayvariance = 0;
    int isgamma = 1;
    float variance = 0.01;
    int minraysamples = 1;
    int maxraysamples = 1;

    float prevlum = .0;
    float var = .0;

    int
    stop_by_variance(const float _lum; const int isample)
    {
	if (!dorayvariance)
	    return 0;

	int i1 = isample + 1;
	if (i1 >= minraysamples)
	    {
		float lum = _lum;
		if (isgamma) lum = sqrt(lum);

		int samplesize;
		float mean;
		float newvar = variance(lum - prevlum,
					mean, samplesize);

		var = (var * isample + newvar) / i1;
		prevlum = lum;

		if (var <= variance*variance)
		    return 1;
	    }
	return 0;
    }
};


struct SamplingFactory
{
    float variance;
    int dorayvariance;
    int isgamma;
    int minraysamples;
    int maxraysamples;

    void init()
    {
	renderstate("object:variance", variance);
	renderstate("object:dorayvariance", dorayvariance);
	renderstate("light:maxraysamples", maxraysamples);
	renderstate("light:minraysamples", minraysamples);
	string colorspace;
	renderstate("renderer:colorspace", colorspace);
	isgamma = colorspace == "gamma";
    }

    VarianceSampler
    getsampler(int quality; int nsamples)
    {
	VarianceSampler ret;

	// Sampling quality
	//	0 - Auto (dorayvariance)
	//	1 - minraysamples
	//	2 - maxraysamples
	//	3 - _nsamples
	if (quality == 0)
	    {
		ret.dorayvariance = 1;
		ret.maxraysamples = maxraysamples;
		ret.minraysamples = minraysamples;
		ret.isgamma = isgamma;
		ret.variance = variance;
	    }
	else
	    {
		ret.dorayvariance = 0;

		if (quality == 1)
		    ret.maxraysamples = minraysamples;
		else if (quality == 2)
		    ret.maxraysamples = maxraysamples;
		else
		    ret.maxraysamples = nsamples;
	    }

	return ret;
    }
};

#endif	//__phy_utils__
