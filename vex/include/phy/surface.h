// This may look like -*- C -*- code, but it is really Houdini Vex
//
//	surface.h - phySurface for Mantra,
//	Physical based, easy to use, compact surface ubershader.
//	This is a part of phyShader for Mantra.
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
#include <phy/sss.h>


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
	      int shadow;
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

    vector sp;
    if (!getsmoothP(sp, -v))
    	sp = p;

    START_ILLUMINANCE;

    vector _tmpDFS = .0;
    vector _tmpSPC = .0;
    vector _tmpTRN = .0;
    vector _tmpSSS = .0;
    float pdfDFS = 0;
    float pdfSPC = 0;
    float pdfTRN = 0;
    float pdfSSS = 0;

    START_SAMPLING("nextpixel");
    SET_SAMPLE;

    SAMPLE_LIGHT(sp, nfN);

    // Diffuse reflection
    if (enableDFS)
	{
	    float weight = 0;
	    vector eval = eval_bsdf(f_DFS, v, l, weight, PBR_DIFFUSE_MASK);
	    _tmpDFS += cl * eval * weight;
	    pdfDFS += weight;
	}

    // Specular reflection
    if (enableSPC)
	{
	    float weight = 0;
	    vector eval = eval_bsdf(f_SPC, v, l, weight, PBR_REFLECT_MASK);
	    _tmpSPC += cl * eval * weight;
	    pdfSPC += weight;
	}
    // Refraction
    if (enableTRN)
	{
	    SAMPLE_LIGHT(pTRN, nbN);

	    float weight = 0;
	    vector eval = eval_bsdf(f_TRN, v, l, weight, PBR_REFRACT_MASK);
	    _tmpTRN += cl * eval * weight;
	    pdfTRN += weight;
	}

    // SSS of thin objects
    if (enableSSS)
    	{
    	    SAMPLE_LIGHT(pSSS, nbN);

	    float weight = 0;
	    vector eval = eval_bsdf(f_SSS, v, l, weight, PBR_DIFFUSE_MASK);
	    _tmpSSS += cl * eval * weight;
	    pdfSSS += weight;
    	}

    END_LOOP;			// SAMPLING

    dfs = factorDFS * _tmpDFS / pdfDFS;
    spc = factorSPC * _tmpSPC / pdfSPC;
    trn = factorTRN * _tmpTRN / pdfTRN;
    sss = factorSSS * _tmpSSS / pdfSSS;

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
			  "ray:length", raylength,	\
			  variable, hitCf,		\
			  "variancevar", hitCf)

#define INIT_TRACE vector hitCf = .0;			\
    vector eval = .0;					\
    int mask = PBR_REFRACT_MASK;			\
    if (raystyle == "reflect") mask = PBR_REFLECT_MASK; \
    string sfilter = oblend ? "opacity" : "closest";	\
    int doabs = (max(absty) > .0);			\
    float raylength = .0

// Delta function case
vector
raytrace(vector p, dir;
	 float maxdist;
	 int oblend, style;
	 string raystyle, scope, variable;
	 vector absty;
	 int dosss;
	 RayMarcher sss_single;
	 export vector sss)
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
			    {
				eval = hitCf * exp(-raylength * absty);

				if (dosss)
				    sss += sss_single->eval(p, dir,
							    raylength);
			    }
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
	 vector absty;
	 int dosss;
	 RayMarcher sss_single;
	 export vector sss)
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
			       "ray:length", raylength,
			       "ray:direction", raydir)
			    {
				tmp += hitCf * exp(-raylength * absty);

				if (dosss)
				    sss += sss_single->eval(p, raydir,
							    raylength);
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
		sss /= samples;
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
#define AVERAGE if (summ > 0) eval /= summ

#define VARIANCEAA if (vsampler->interrupt(max(eval)/summ, _i)) break;

#define FINALIZE_SAMPLING CONTRIBUTE; VARIANCEAA; END_LOOP; AVERAGE

vector
raytrace(bsdf f;
	 vector p, v;
	 float maxdist;
	 int sid, oblend, style;
	 string raystyle, scope, variable;
	 vector absty;
	 VarianceSampler vsampler;
	 int dosss;
	 RayMarcher sss_single;
	 export vector sss)
{
    INIT_TRACE;

    int type = 0;
    vector dir = 0;
    vector brdf = 0;
    vector tmp = 0;
    float pdf = 0;
    float summ = .0;
    int samples = vsampler.maxraysamples;

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
			if(trace(p, dir, Time,
				 "scope", scope,
				 "samplefilter", "closest",
				 "maxdist", maxdist,
				 "raystyle", raystyle,
				 "ray:length", raylength,
				 variable, hitCf,
				 "variancevar", hitCf))
			    {
				vector absrp = exp(-raylength * absty);

				if (dosss)
				    sss += pdf * sss_single->eval(p, dir,
								  raylength);

				tmp = max(hitCf, .0) * absrp;
			    }
		    }
		else
		    {
			if (TRACE_GATHER)
			    tmp = hitCf;
			else
			    tmp = resolvemissedray(dir, Time, mask);
		    }

		if (summ > 0)
		    sss /= summ;
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
    vector eval = .0;				\
    sss = .0

// volume absorption
// Delta case
vector
absorption(vector p, dir, kabs;
	   float maxdist;
	   string scope;
	   int dosss;
	   RayMarcher sss_single;
	   export vector sss)
{
    INIT_ABSRP;

    if(trace(p, dir, Time,
	     TRACE_FLAGS(1),
	     "ray:length", raylength))
	{
	    eval = exp(-raylength * kabs);

	    if (dosss)
		sss = sss_single->eval(p, dir,
				      raylength);
	}

    return eval;
}

// Cone case
vector
absorption(vector p, dir, kabs;
	   float maxdist, angle;
	   int samples;
	   string scope;
	   int dosss;
	   RayMarcher sss_single;
	   export vector sss)
{
    INIT_ABSRP;

    gather(p, dir,
	   TRACE_FLAGS(samples),
	   "angle", angle,
	   "ray:length", raylength)
	{
	    eval += exp(-raylength * kabs);

	    if (dosss)
		sss = sss_single->eval(p, dir,
				      raylength);
	}

    sss /= samples;
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


// optimize X absorption
#define OPTABS(X) if (X) { rtAbsty = _absty; rtSSS = allowsinglesss ?  _sssca : .0; }


// The main surface routine
void
physurface(int conductor;
	   int thin;
	   float thickness;	// Thickness of thin sheet if thin
	   float alb;		// Albedo
	   vector clrsurf;	// Surface color
	   vector2 iort;	// Pair of refraction indicies of transmission medium
	   float iori;		// refraction index of incident medium
	   int enableDFS, enableSPC, enableTRN, enableSSS;
	   float weightDFS, weightSPC, weightTRN, weightSSS;
	   float roughDFS; 	// Oren-Nayar roughness
	   float roughSPC;	// GTR alpha parameter
	   float gamma;		// GTR gamma parameter
	   float anisobias;	// Anisotropy reflection/refraction bias
	   vector sssca;	// SSS albedo
	   vector sssdf;	// SSS diffuse
	   vector absty;	// Absorption coefficient
	   float g;		// Scattring phase
	   int dispersion;
	   int styleSPC, styleTRN;   // How to perform reflection/refraction
	   int oblendSPC, oblendTRN; // Opacity blending
	   int tsquality;	// Ray-tracing sampling quality
	   int msquality; 	// Multiple scattering sampling quality
	   int vsquality;	// Single scattering sampling quality
	   int tsamples;	// Number of ray-tracing samples
	   int vsamples;	// Number of single scattering samples
	   int msamples;	// Number of multiple scattering samples
	   int mdisablesecondary; // Disable secondary multiple scattering
	   int vdisablesecondary; // Disable secondary single scattering
	   int shadow;		// Receive shadows
	   int empty;
	   int useF;		// Use BSDF to compute reflection/refraction
	   float depthimp;	// Depth importance
	   float maxdist;	// Maximum tracing distance
	   vector p;		// Point position
	   vector n;		// Point normal
	   vector ii;		// Incident ray from eye to point
	   vector tangent;	// Surface derivative
	   vector sellmeierB, sellmeierC; // Sellmeier's coefficients
	   string gvarSPC, gvarTRN;	  // Gather variables
	   float density;		  // Single scattering density
	   string sscope;		  // Multiple scattering object scope
	   float curvature;		  // SSS sampling parameter
	   string lightmasksss;		  // SSS lightmask
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
    int sid = renderengine == "micropoly" ? newsampler() : newsampler(SID);


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
    float etat = max(1.0000016, iort.x);

    // Dispersion
    vector dtint = 1.;
    if (dispersion)
	{
	    float wl = samplewl(rand(sid));
	    etat = sellmeier(wl, sellmeierB, sellmeierC);
	    dtint = wl2rgb(wl);
	}

    float
	etai = iori,
	etak = iort.y,
	eta = etai/etat;

    // Geometry atenuation factor
    float
	gafmask = 1.,
	gafrefr = 1.;

    // Treat surface as perfect smooth
    // if specular rougness is below this value
    int smooth = roughSPC <= SMOOTH_THRESHOLD;
    
    // Ray-tracing spread angle
    float angle = .0;

    // Absorption coefficient, scalar
    float kabs = max(absty);

    int doAbs = kabs > .0;

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

    // Normalized SSS color 
    vector _sssca = sssca / ALONE_VEC(sssca);
    // Scattering coefficint
    vector _sca = density * invert_hue(sssca);

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

    // Is the total internal reflection case
    int internal = rdir == tdir;

    int solid = !thin;
    int thick = thin && thickness > .0;


    int depth = max(getraylevel(), getglobalraylevel());
    SamplingFactory sfactory;
    sfactory->init(depth, depthimp);

    // raytracing variance aa
    VarianceSampler tvsampler = sfactory->getsampler(tsquality, tsamples);
    // multiple scattering variance aa
    VarianceSampler mvsampler = sfactory->getsampler(msquality, msamples);

    // Ray-tracing scope
    string scopeSPC = "scope:default";
    string scopeTRN = "scope:default";

    bsdf f_DFS = cvex_bsdf("phy_diffuse_eval",
    			   "phy_diffuse_sample",
    			   "label", "diffuse",
    			   "n", nfN,
    			   "sigma", roughDFS,
    			   "eta", enableSPC ? eta : 1.);
    bsdf f_SPC = specular(rdir);
    bsdf f_TRN = specular(tdir);
    bsdf f_SSS = cvex_bsdf("phy_diffuse_eval",
			   "phy_diffuse_sample",
			   "label", "diffuse",
			   "n", nbN);
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
	msamples	&&
	!(mdisablesecondary && depth);

    int allowsinglesss =
	allowSSS	&&
	vsamples	&&
    	!(vdisablesecondary && depth);

    int translucent =
	enableSSS	&&
	thin;

    // Single scattering
    RayMarcher sss_single;
    sss_single.ca = _absty;
    sss_single.f = f_VOL;
    sss_single.sid = sid;
    sss_single.depth = depth;
    sss_single.depthimp = depthimp;
    sss_single.doshadow = shadow;
    sss_single.lightmask = lightmasksss;
    sss_single.vsampler = sfactory->getsampler(vsquality, vsamples);

    // disable separate absorption and single scattering
    // for Raytrace/Micropoly renderers
    int inside = enter || internal;

    int isRTMP = (renderengine == "micropoly" ||
		  renderengine == "raytrace");

    int oAbs = doAbs && isRTMP;

    // refraction or internal
    int oAbsTRN = oAbs && inside;

    // internal reflection will perform, drop specular
    int oAbsSPC = oAbs && !inside;

    vector rtAbsty = .0;
    vector rtSSS = .0;


    // Get tracing masks for non-gather tracing
    if(!renderstate("object:reflectmask", scopeSPC))
	scopeSPC = "scope:default";
    if(!renderstate("object:refracttmask", scopeTRN))
	scopeTRN = "scope:default";

    if (empty && solid)
	{
	    if (inside)
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


    // BSDFs of reflections and refractions
    if(smooth)
	{
	    f_SPC = specular(rdir);
	    f_TRN = specular(tdir);
	}
    else
	{
	    float sigma = max(roughSPC*roughSPC, SMOOTH_THRESHOLD);

	    angle = atan(sigma);

	    vector2 alpha = anisorough(sigma, anisobias);
	    float beta = anisorough_i(v, tU, tV, alpha);

	    gafmask = gaf(dot(v, nfN), beta);
	    gafrefr = gaf(dot(tdir, ni), beta);

	    f_SPC = get_ggr("reflect",
			    nfN, tU, tV,
			    alpha, gamma, 1.);

	    f_TRN = get_ggr("refract",
			    n, tU, tV,
			    alpha, gamma, thin ? 1. : eta);
	}

    // Ray-tracing specular
    if (allowSPC && styleSPC && isRTMP)
	{
	    OPTABS(oAbsSPC);

	    if (smooth)
		traceSPC = raytrace(p, rdir,
				    maxdist,
				    oblendSPC, styleSPC,
				    "reflect", scopeSPC, gvarSPC,
				    rtAbsty,
				    allowsinglesss,
				    sss_single,
				    singlescattering);
	    else if (useF)
		traceSPC = raytrace(f_SPC,
				    p, v,
				    maxdist,
				    sid, oblendSPC, styleSPC,
				    "reflect", scopeSPC, gvarSPC,
				    rtAbsty,
				    tvsampler,
				    allowsinglesss,
				    sss_single,
				    singlescattering);
	    else
		traceSPC = raytrace(p, rdir,
				    angle, maxdist,
				    tvsampler.maxraysamples,
				    oblendSPC, styleSPC,
				    "reflect", scopeSPC, gvarSPC,
				    rtAbsty,
				    allowsinglesss,
				    sss_single,
				    singlescattering);
	}

    // Refraction
    if (allowTRN && styleTRN)
	{
	    int do_trace = styleTRN && isRTMP;

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
	    else if (doAbs && !oAbs)
		{
		    // Absorption and single scattering for PBR
		    vector abstmp = 1.;
		    vector tmpsss = .0;

		    if (inside)
			absdir = tdir;
		    else if (enableSPC)
			absdir = rdir;

		    if (smooth)
		    	abstmp = absorption(p, absdir,
		    			    _absty,
		    			    maxdist,
		    			    scopeTRN,
					    allowsinglesss,
		    			    sss_single,
		    			    tmpsss);
		    else
		    	abstmp = absorption(p, absdir,
		    			    _absty,
		    			    maxdist, angle,
		    			    tvsampler.maxraysamples,
		    			    scopeTRN,
					    allowsinglesss,
		    			    sss_single,
		    			    tmpsss);

		    singlescattering = tmpsss;

		    if (inside)
			absTRN = abstmp;
		    else if(enableSPC)
			absSPC = abstmp;
		}

	    if (do_trace)
		{
		    OPTABS(oAbsTRN);

		    if (smooth)
			traceTRN = raytrace(pTRN, tdir,
					    maxdist,
					    oblendTRN, styleTRN,
					    "refract", scopeTRN, gvarTRN,
					    rtAbsty,
					    allowsinglesss,
					    sss_single,
					    singlescattering);
		    else if (useF)
			traceTRN = raytrace(f_TRN,
					    pTRN, v,
					    maxdist,
					    sid, oblendTRN, styleTRN,
					    "refract", scopeTRN, gvarTRN,
					    rtAbsty,
					    tvsampler,
					    allowsinglesss,
					    sss_single,
					    singlescattering);
		    else
			traceTRN = raytrace(pTRN, tdir,
					    angle, maxdist,
					    tvsampler.maxraysamples,
					    oblendTRN, styleTRN,
					    "refract", scopeTRN, gvarTRN,
					    rtAbsty,
					    allowsinglesss,
					    sss_single,
					    singlescattering);
		}
	}


    // The base factors
    factorDFS = clrDFS * kDFS;
    factorSPC = clrSPC * kSPC * fr * gafmask;
    factorTRN = kTRN * ft * gafrefr * dtint;
    factorSSS = kSSS * sssdf;

    // Translucency
    if (thin)
	{
	    factorSSS *= _sssca;
	    factorTRN *= clrTRN * _absTRN;
	}

    // Lighting
    illum_surface(p, pTRN, p + nbN * thickness,
    		  nfN, nbN,
    		  v,
    		  thick,
    		  thickness,
    		  depthimp,
		  shadow,
    		  sid,
    		  depth,
    		  enableDFS, allowSPC, allowTRN, translucent,
    		  f_DFS, f_SPC, f_TRN, f_SSS,
    		  factorDFS, factorSPC, factorTRN, factorSSS,
    		  fullDFS, fullSPC, fullTRN, fullSSS);

    // Compute multiple scattering
    if (allowmultisss)
    	fullSSS = sss_multi(p, n,
			    sssca,
			    eta,
			    sid,
			    sscope,
			    shadow,
			    curvature,
			    lightmasksss,
			    depth, depthimp,
			    mvsampler)
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
		sh *= shadowabs(p, allowsinglesss ? _sca : _absty);

	    opacity = 1. - sh;
	}

    beauty = fullDFS + fullSPC + fullTRN + fullSSS;
    all = array(fullDFS, fullSPC, fullTRN, fullSSS);
    f = f_DFS + f_SPC + f_TRN + f_SSS;
}


#endif // __physurface__
