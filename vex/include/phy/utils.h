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

#define START_SAMPLING						\
    float sx, sy;						\
    for (int idx = 0; idx < samples; ++idx) {			\
    nextsample(sid, sx, sy, "mode", "nextpixel");

#define SET_SAMPLE	vector sample = set(sx, sy, .0);

#define END_LOOP	}


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


#endif	//__phy_utils__
