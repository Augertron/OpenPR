//////////////////////////////////////////////////////////////////////////////////////////////////
// Author:		Styx(Zhenhui Xu)
// Version:		0.1
// Date:		May 20, 2009 
// Description: Heap sort algorithms.         
// Others:      Defination of the class Image.
// 
// Copyright(C) 2009 OpenPR
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, 
// are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright notice, this
//       list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright notice, 
//       this list of conditions and the following disclaimer in the documentation 
//       and/or other materials provided with the distribution.
//     * Neither the name of the OpenPR nor the names of its contributors may
//       be used to endorse or promote products derived from this software without 
//       specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT 
// SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED 
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR 
// BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN 
// ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//////////////////////////////////////////////////////////////////////////////////////////////////


#include "Image.h"


Image::Image(void)
{
	m_pixels = 0;
}

Image::Image(int xDim, int yDim)
{
	this->Allocate(xDim, yDim);
}

Image::Image(const Image &other)
{
	this->Allocate(other.m_xDim, other.m_yDim);

	for (int y=0; y<m_yDim; y++)
	{
		for (int x=0; x<m_xDim; x++)
		{
			m_pixels[y][x] = other(x, y);
		}
	}
}

Image::~Image(void)
{
	DeAllocate();
}

Image& Image::operator= (const Image& other)
{
	if (this == &other)
	{
		return *this;
	}

	if ( (m_xDim != other.GetXDim()) || (m_yDim != other.GetYDim()) || m_pixels == 0 )
		ReAllocate(other.GetXDim(), other.GetYDim());

	for (int y=0; y<m_yDim; y++)
	{
		for (int x=0; x<m_xDim; x++)
		{
			// m_pixels[y][x] = other.m_pixels[y][x];
			m_pixels[y][x] = other(x,y);
		}
	}
	return *this;
}


void Image::ReAllocate(int xDim, int yDim)
{
	this->DeAllocate();
	this->Allocate(xDim, yDim);
}

void Image::Allocate(int xDim, int yDim)
{
	assert((xDim > 0) && (yDim > 0));

	m_xDim = xDim;
	m_yDim = yDim;

	m_pixels = new double *[m_yDim];
	if (m_pixels == NULL)
	{
		FatalError("Image -- Allocating memory fails!");
	}
	for (int y=0; y<m_yDim; y++)
	{
		m_pixels[y] = new double[m_xDim];
		if (m_pixels[y] == NULL)
		{
			FatalError("Image -- Allocating memory fails!");
		}
	}
}

void Image::DeAllocate(void)
{
	if (m_pixels != NULL)
	{
		for (int y=0; y<m_yDim; y++)
		{
			delete [] m_pixels[y];
		}
		delete []m_pixels;
	}
}

/********************************************************************
   So we can just use image(x, y) to indicate the pixel value 
   at position (x, y).
********************************************************************/
double & Image::operator ()(int x, int y)
{
	assert((x >= 0) && (x < m_xDim) && (y >= 0) && (y < m_yDim));

	return m_pixels[y][x];
}

const double& Image::operator() (int x, int y) const
{
	assert((x >= 0) && (x < m_xDim) && (y >= 0) && (y < m_yDim));

	return m_pixels[y][x];
}

int Image::GetXDim(void) const
{
	return m_xDim;
}

int Image::GetYDim(void) const
{
	return m_yDim;
}

/********************************************************************
   Find the minimum to maximum range, then stretch and limit those 
   to exactly 0.0 to 1.0. 
   If the minimum and maximum values are equal, no normalization 
   takes place.
********************************************************************/
void Image::Normalize()
{
	double min = 1.0;
	double max = 0.0;
	int x, y;

	for (y=0; y<m_yDim; ++y)
	{
		for (x=0; x<m_xDim; ++x)
		{
			if (min > m_pixels[y][x])
			{
				min = m_pixels[y][x];
			}
			
			if (max < m_pixels[y][x])
			{
				max = m_pixels[y][x];
			}
		}
	}
	
	if (min == max)
	{
		return;
	}

	double diff = max - min;

	for (y=0; y<m_yDim; ++y)
	{
		for (x=0 ; x<m_xDim; ++x)
		{
			m_pixels[y][x] = (m_pixels[y][x] - min) / diff;
		}
	}
}

/********************************************************************
   Downscale the image from size (x, y) to (x / 2, y / 2), 
   sampling pixels by a factor of 2.
********************************************************************/
Image Image::HalfScale(void)
{
	if ((m_xDim / 2 == 0) || (m_yDim / 2 == 0))
	{
		//return *this;
		FatalError("Too small image size, cannot half.");
	}

	Image res(m_xDim / 2, m_yDim / 2);
	
	for (int y=0; y<res.m_yDim; y++)
	{
		for (int x=0; x<res.m_xDim; x++)
		{
			res.m_pixels[y][x] = this->m_pixels[2 * y][2 * x];
		}
	}

	return res;
}

/********************************************************************
   Double the image from size (x, y) to (x * 2 - 1, y * 2 - 1), 
   using linear interpolation.
********************************************************************/
Image Image::DoubleScale(void)
{
	if ((m_xDim <= 1) || (m_yDim <= 1))
	{
		//return *this;
		FatalError("Too small image size for doubling.");
	}

	Image res(m_xDim * 2 - 1, m_yDim * 2 - 1);

	int x, y;
	for (y=0; y<(m_yDim-1); y++)
	{
		for (x=0; x<(m_xDim-1); x++)
		{
			res(x * 2, y * 2) = m_pixels[y][x];
			res(x * 2 + 1, y * 2) = 0.5 * 
				(m_pixels[y][x] + m_pixels[y][x + 1]);
			res(x * 2, y * 2 + 1) = 0.5 *
				(m_pixels[y][x] + m_pixels[y + 1][x]);
			res(x * 2 + 1, y * 2 + 1) = 0.25 *
				(m_pixels[y][x] + m_pixels[y + 1][x + 1] +
				m_pixels[y][x + 1] + m_pixels[y + 1][x]);
		}
	}

	for (y=0; y<(m_yDim-1); y++)
	{
		res(m_xDim * 2 - 2, y * 2) = m_pixels[y][m_xDim - 1];
		res(m_xDim * 2 - 2, y * 2 + 1) = 0.5 * 
			(m_pixels[y][m_xDim - 1] + 
			m_pixels[y + 1][m_xDim - 1]);
	}

	for (int x=0; x<(m_xDim-1); x++)
	{
		res(x * 2, m_yDim * 2 - 2) = m_pixels[m_yDim - 1][x];
		res(x * 2 + 1, m_yDim * 2 - 2) = 0.5 * 
			(m_pixels[m_yDim - 1][x] + 
			m_pixels[m_yDim - 1][x + 1]);
	}

	res(m_xDim * 2 - 2, m_yDim * 2 - 2) = 
		m_pixels[m_yDim - 1][m_xDim - 1];

	return res;
}

/*
Image operator -(Image &img1, Image &img2)
{
	if ((img1.m_xDim != img2.m_xDim) || (img1.m_yDim != img2.m_yDim))
	{
		cerr << "Images have different sizes, can not subtract.\n";
		exit(1);
	}

	Image res(img1.m_xDim, img1.m_yDim);

	for (int y=0; y<img1.m_yDim; y++)
	{
		for (int x=0; x<img1.m_xDim; x++)
		{
			res(x, y) = img1(x, y) - img2(x, y);
		}
	}

	return res;
}
*/

Image operator- (Image& img1, Image& img2)
{
	int dimX = img1.GetXDim();
	int dimY = img1.GetYDim();

	if ((dimX != img2.GetXDim()) || (dimY != img2.GetYDim()))
	{
		FatalError("Images have different sizes, can not subtract.");
	}

	Image res(dimX, dimY);

	for (int y=0; y<dimY; y++)
	{
		for (int x=0; x<dimX; x++)
		{
			res(x, y) = img1(x, y) - img2(x, y);
		}
	}

	return res;
}
