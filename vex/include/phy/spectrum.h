// This may look like -*- C -*- code, but it is really Houdini Vex
// 
//	spectrum.h - Wavelength depended routines.
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


#ifndef __phy_spectrum__
#define __phy_spectrum__


#include <phy/utils.h>


// All wavelength are measured in micrometers,
// except specified ones


#define LOWER	0.400		// The lower limit of visible spectrum
#define HIGHER	0.680		// The higher limit of visible spectrum
#define RANGE (HIGHER - LOWER)

#define CIE_REC_709	set(3.2404542, -1.5371385, -0.4985314,	\
			    -0.9692660, 1.8760108, 0.0415560,	\
			    0.0556434, -0.2040259, 1.0572252);

// White for the given range
#define WHITE set(0.68522342, 0.29870292, 0.42072795)

// Wavelength to CIE XYZ
// 
// Based on:
// Chris Wyman, Peter-Pike Sloan, Peter Shirley`s
// Simple Analytic Approximations to the CIE XYZ Color Matching Functions
vector
wl2xyz(float iwl)
{
    // From micrometers to nanometers
    float wl = iwl * 1000.;

    float x1 = (wl-442.0) * ((wl<442.0) ? 0.0624 : 0.0374);
    float x2 = (wl-599.8) * ((wl<599.8) ? 0.0264 : 0.0323);
    float x3 = (wl-501.1) * ((wl<501.1) ? 0.0490 : 0.0382);
    float x
	= 0.362 * exp(-0.5*x1*x1)
	+ 1.056 * exp(-0.5*x2*x2)
	- 0.065 * exp(-0.5*x3*x3);

    float y1 = (wl-568.8) * ((wl<568.8) ? 0.0213 : 0.0247);
    float y2 = (wl-530.9) * ((wl<530.9) ? 0.0613 : 0.0322);
    float y
	= 0.821 * exp(-0.5*y1*y1)
	+ 0.286 * exp(-0.5*y2*y2);

    float z1 = (wl-437.0) * ((wl<437.0) ? 0.0845 : 0.0278);
    float z2 = (wl-459.0) * ((wl<459.0) ? 0.0385 : 0.0725);
    float z
	= 1.217 * exp(-0.5*z1*z1)
	+ 0.681 * exp(-0.5*z2*z2);

    return set(x, y, z);
}


// CIE XYZ to RGB
vector
xyz2rgb(vector xyz)
{
    vector eval = xyz * CIE_REC_709;
    return FIT01(eval);
}


// Wavelength to RGB
vector
wl2rgb(float wl)
{
    return xyz2rgb(wl2xyz(wl)) / WHITE;
}


// Get IOR for choosen wavelength
// by Sellmeier equation
float
sellmeier(float wl;
	  vector Bk, Ck)
{
    float wl2= wl * wl;
    return sqrt(1.
                + Bk.x * wl2 / (wl2 - Ck.x)
                + Bk.y * wl2 / (wl2 - Ck.y)
                + Bk.z * wl2 / (wl2 - Ck.z));
}


// Sample wavelength
float
samplewl(float sx)
{
    return LOWER + RANGE * sx;
}


#endif
