//////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2009 OpenPR
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of OpenPR nor the names of its 
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL HOLDER AND CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
////////////////////////////////////////////////////////////////////////////

// Image.h -- Declaration of the class Image.

#ifndef IMAGE_H
#define IMAGE_H

/*
#include <cstdlib>
#include <iostream>
#include <cassert>

using namespace std;
*/
#include "SiftUtil.h"


class Image
{
public:
	Image();
	Image(int xDim, int yDim);
	Image(const Image &other);
	~Image();
	Image & operator =(const Image &other);

	void ReAllocate(int xDim, int yDim);
	void Allocate(int xDim, int yDim);
	void DeAllocate();
	
	double& operator() (int x, int y);
	const double& operator() (int x, int y) const;

	int GetXDim(void) const;
	int GetYDim(void) const;

	void Normalize();
	Image HalfScale();
	Image DoubleScale();
	//friend Image operator -(Image &img1, Image &img2);

private:
	int m_xDim; // width of the image.
	int m_yDim; // height of the image.
	double **m_pixels; // pixel values of the image.
};


Image operator -(Image &img1, Image &img2);


#endif
