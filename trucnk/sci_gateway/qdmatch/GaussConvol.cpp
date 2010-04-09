//////////////////////////////////////////////////////////////////////////////////////////////////
// Author:		Styx(Zhenhui Xu)
// Version:		0.1
// Date:		May 20, 2009 
// Description: This file is used to compute Gaussian images.         
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


#include "GaussConvol.h"


/********************************************************************
   Compute Gaussian convolution of the image.
********************************************************************/

// Modified by styx. (Dec 12, 2007)
Image GaussConvol(const Image& src, double sigma)
{
	Image dst(src.GetXDim(), src.GetYDim());

   // Size of the Gaussian convolution kernel is 
   // (1 + 2 * (int)(3.0 * sigma + 0.5)).
	int kernelSize = 1 + 2 * static_cast<int>(3.0 * sigma + 0.5);

	double *kernel = new double[kernelSize];
	
	GaussKernel(kernelSize, sigma, kernel);

	if (kernel != NULL)
	{
		Convolve2D(dst, src, kernelSize, kernel);
		delete []kernel;
	}
	else
	{
		// dst = src;
		return src;
	}
	return dst;
}

/********************************************************************
 *             Calculate 1-D Gaussian Kernel.
 ********************************************************************/
void GaussKernel(int kernelSize, double sigma, double* kernel)
{
	//kernelSize = 1 + 2 * (int)(3.0 * sigma + 0.5);
	int kernelMid = kernelSize / 2;

	//double *kernel = new double[kernelSize];

	double sigma2sq = 1.0 / (2.0 * sigma * sigma);
	double normFactor = 1.0 / (sqrt(2.0 * PI) * sigma);

	for (int n=0; n<kernelSize; n++)
	{
		double expFactor = - ((n-kernelMid) * (n-kernelMid)) * sigma2sq;

		kernel[n] = exp(expFactor) * normFactor;
	}
}

/********************************************************************
   Compute the convolution of the image with the convolution kernel
   in both vertical and horizontal direction.
********************************************************************/
void Convolve2D(Image& dst, const Image& src, 
				int kernelSize, const double* kernel)
{
	Image tmp(src.GetXDim(), src.GetYDim());

	Convolve1D(tmp, src, kernelSize, kernel, VERTICAL);
	Convolve1D(dst, tmp, kernelSize, kernel, HORIZONTAL);
}

/********************************************************************
   Separable convolution, compute in one direction.
********************************************************************/
void Convolve1D(Image& dst, const Image& src, 
				int kernelSize, const double *kernel, DIRECTION dir)
{
	int maxN; // outer loop maximum index.
	int maxP; // inner loop maximum index.

	if (dir == VERTICAL)
	{
		maxN = src.GetXDim();
		maxP = src.GetYDim();
	}
	else if (dir == HORIZONTAL)
	{
		maxN = src.GetYDim();
		maxP = src.GetXDim();
	}
	else
	{
		FatalError("Convolve1D -- Invalid DIRECTION!");
	}

	int idx = (kernelSize - 1) / 2;

	/////
	//clock_t start = clock();

	for (int n=0; n<maxN; n++)
	{
		for (int p=0; p<maxP; p++)
		{
			int k = p + idx;
			//int jMin = Max(0, k + 1 - kernelSize);
			//int jMax = Min(k + 1, maxP);
			int jMin = Maximun(0, k + 1 - kernelSize);
			int jMax = Minimun(k + 1, maxP);
			double val = 0.0;

			if (dir == VERTICAL)
			{
				for (int j=jMin; j<jMax; j++)
				{
					val += src(n, j) * kernel[k - j];
				}

				dst(n, p) = val;
			}
			else
			{
				for (int j=jMin; j<jMax; j++)
				{
					val += src(j, n) * kernel[k - j];
				}

				dst(p, n) = val;
			}
		}
	} // end of for.

	/////
	//clock_t finish = clock();
	//double time = (double)(finish - start) / CLOCKS_PER_SEC;
	//cout << "for uses " << time << " s." << endl;
}
