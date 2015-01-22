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


#include <phy/sss.h>
#include <phy/utils.h>
#include <phy/shading.h>
#include <phy/spectrum.h>
#include <phy/microfacet.h>


#define MAX_ROUGH	0.3


// Return number of samples
// by mode.
// mode choices:
//	0 - Auto
//	1 - Minimum
//	2 - Maximum
//	3 - Manual
// If mode = Auto use sigma to compute mean
int
sampling_quality(int mode;
		 int number;
		 float sigma)
{
    int
	eval = number,
	maxsamples = 1,
	minsamples = 1;

    renderstate("light:maxraysamples", maxsamples);
    renderstate("light:minraysamples", minsamples);

    if (mode != 3)
	if (mode == 0)
	    eval = floor(lerp(minsamples, maxsamples,
			      max(sigma / MAX_ROUGH, MAX_ROUGH)));
	else if (mode == 1)
	    eval = minsamples;
	else
	    eval = maxsamples;

    return eval;
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
	   int styleSPC, styleTRN;   // How to perform reflection/refraction (see phy_shading.h)
	   int oblendSPC, oblendTRN; // Opacity blending
	   int squality;	// Sampling quality (see sampling_quality)
	   int tsamples;	// Number of ray-tracing samples
	   int vsamples;	// Number of single scattering samples
	   int dsamples;	// Number of dispersion samples
	   int ssamples;	// Number of multiple scattering samples
	   int empty;
	   int accurateabs;	// Accurate absorption
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
    float
	etat = max(1.0000016, ior.x),
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
    int issmooth = roughSPC <= SMOOTH_THRESHOLD;

    // When reflections/refraction are anisotropic
    int isanisotropic = (anisobias != .0) && !issmooth;

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

    // Dispersions
    vector dcolors[], tdirs[];

    // Ray is directed inside object
    int	enter = dot(n, ni) < .0;

    int depth = getraylevel() + getglobalraylevel();

    // Copied variables starts with "_"
    int _tsamples = sampling_quality(squality, tsamples, sigma);
    int _vsamples = vsamples;
    int _ssamples = ssamples;
    int _dsamples = dsamples;

    // Is the total internal reflection case
    int internal = rdir == tdir;

    int solid = !thin;
    int thick = thin && thickness > .0;

    // Recompute number of samples by the depth importance
    if (depth && depthimp < 1.)
	{
	    float factor = pow(depthimp, depth);

	    _tsamples = FLOOR_ALONE(tsamples * factor);
	    _dsamples = FLOOR_ALONE(dsamples * factor);
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

    string renderengine;
    renderstate("renderer:renderengine", renderengine);

    int sid = renderengine == "micropoly" ? newsampler() : SID;

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
	
    // When dispersion is allowed
    int allowdispersion =
	allowTRN	&&
	solid		&&
	(_dsamples > 1)	&&
	issmooth;

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

    if (allowdispersion)
	{
	    dispersions(ni, n,
			sellmeierB, sellmeierC,
			enter,
			sid, _dsamples,
			etai,
			dcolors, tdirs);

	    // For Fresnel refletance use
	    // Mid wavelength eta of interior
	    etat = sellmeier(RANGE *.5 + LOWER,
			     sellmeierB, sellmeierC);
	    // Update eta
	    eta = etai/etat;
	}

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
	    fr = cfresnel(ni, nfN, etat, etak);
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
    if(issmooth)
	{
	    f_SPC = specular(rdir);

	    if (allowdispersion)
		{
		    for (int idx = 0; idx < _dsamples; ++idx)
			f_TRN += dcolors[idx] *
			    specular(tdirs[idx]);

		    f_TRN *= 1. / (float)_dsamples;
		}
	    else
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
	if (issmooth)
	    traceSPC = raytrace(p, rdir,
				maxdist,
				oblendSPC, styleSPC,
				"reflect", scopeSPC, gvarSPC);
	else
	    traceSPC = raytrace(p, rdir,
				angle, maxdist,
				_tsamples, oblendSPC, styleSPC,
				"reflect", scopeSPC, gvarSPC);

    // Refraction
    if (allowTRN && styleTRN)
	{
	    if (thin)
		if (thick)
		    thinP(p, ni, nbN, tdir,
			  kabs, eta, thickness,
			  pTRN, _absTRN);
		else
		    {
			if (issmooth)
			    {
				opacity = 1. - kTRN * ft * clrTRN;

				// Don't pathtrace refraction
				// because already used opacity
				f_TRN *= .0;
			    }
		    }
	    else if (kabs > .0)
		{
		    vector abstmp = 1.;
		    vector tmpsss = .0;

		    if (accurateabs && !issmooth)
			{
			    abstmp = absorption(p, absdir,
						_absty,
						maxdist, angle,
						_tsamples,
						scopeTRN);
			}
		    else
			{
			    abstmp = absorption(p, absdir,
						_absty,
						maxdist,
						scopeTRN);
			}


		    if (allowsinglesss)
			tmpsss = raymarch(p, absdir,
					  _absty, clrSSS,
					  g,
					  _vsamples, sid,
					  scopeTRN,
					  depth, depthimp);

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

	    // Setup dispersion
	    if (allowdispersion)
		{
		    for (int idx = 0; idx < _dsamples; ++idx)
			traceTRN += dcolors[idx] *
			    raytrace(pTRN, tdirs[idx],
				     maxdist,
				     oblendTRN, styleTRN,
				     "refract", scopeTRN, gvarTRN);

		    traceTRN /= _dsamples;
		}
	    else
		if (issmooth)
		    {
			// In other cases opacity is used
			if (thick || solid)
			    traceTRN = raytrace(pTRN, tdir,
						maxdist,
						oblendTRN, styleTRN,
						"refract", scopeTRN, gvarTRN);
		    }
		else
		    traceTRN = raytrace(pTRN, tdir,
					angle, maxdist,
					_tsamples, oblendTRN, styleTRN,
					"refract", scopeTRN, gvarTRN);
	}


    // The base factors
    factorDFS = clrDFS * kDFS;
    factorSPC = clrSPC * kSPC * fr * gafmask;
    factorTRN = kTRN * ft * gafrefr;
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
