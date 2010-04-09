//////////////////////////////////////////////////////////////////////////////////////////////////
// Author:		Styx(Zhenhui Xu)
// Version:		0.1
// Date:		May 20, 2009 
// Description: SIFT routine.
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


#include "Sift.h"

#define DEBUG

/********************************************************************
   SIFT function.
   Input an image of PGM format, output a list of keypoints.
   The input parameter 'contrastThresh' can adjust the number of 
   extracted keypoints, usually using 0.01-0.03.

   This algorithm follows the introduction in David Lowe's paper 
   "Distinctive Image Features from Scale-Invariant Keypoints" 
   published in IJCV 2004. The parameters used in this routine 
   are recommended in the paper.

   parameters for building scale space: 
   'preSigma' -- preblur of the original image.
   'startScale' -- subsample of the first octave.
   'scales' -- number of scales sampled per octave. 
               level interval = pow(2.0, 1/scales).
   'octaveSigma' -- prior image smoothing for the initial level 
                    of each octave.
   'minSize' -- size of the images in the last octave.

   parameters for detecting extrema:
   'border' -- the extrema are detected inside the border.

   parameters for accurating keypoints: 
   'maxSteps' -- the largest times of adjustment loop in the step of 
                 accurating keypoints localization.
   'contrastThresh' -- the threshold to discard keypoints with low 
                       contrast. the lower, the more keypoints.
   'edgeRatio' -- the threshold to discard keypoints along the edge.

   parameters for assigning orientations: 
   'binsNum' -- bins number of the orientation histogram.
   'peakThresh' -- threshold with which we find other local peaks.

   parameters for computing descriptors: 
   'gridDim' -- grid dimension of the descriptor.
   'dirNum' -- number of discretized direction.
   'gridSpace' -- grid spacing of the descriptor.
   'illuThresh' -- threshold to avoid illumination change.
********************************************************************/


vector<Keypoint> SIFT(Image &srcImg, double contrastThresh)
{
	cout << "Finding keypoints..." << endl;

	/////
	//clock_t start, finish;
	//double time;

	/////
	//start = clock();

	// normalize input image.
	srcImg.Normalize();
	//WritePGMFile("norm.pgm", srcImg);

	double startScale;
	int maxImgSize = 800;
	int max2Size = 2 * maxImgSize;

	Image preImage;

	if (srcImg.GetXDim() <= maxImgSize || srcImg.GetYDim() <= maxImgSize)
	{
#ifdef DEBUG
		cout << "double image size." << endl;
#endif
		preImage = srcImg.DoubleScale();
		startScale = 0.5; // original is 1.0, after doubling, is 0.5.
	}
	else if (srcImg.GetXDim() > max2Size || srcImg.GetYDim() > max2Size)
	{
		// FatalError("Input image is too large, please resize first.");
		cout << "no doubling." << endl;
		preImage = srcImg;
		startScale = 1.0;
	}
	else
	{
#ifdef DEBUG
		cout << "no doubling." << endl;
#endif
		preImage = srcImg;
		startScale = 1.0;
	}

	/////
	//finish = clock();
	//time = (double)(finish - start) / CLOCKS_PER_SEC;
	//cout << "double image uses " << time << " s." << endl;

	double preSigma = 1.6;

	//double startScale = 0.5; // original is 1.0, after doubling, is 0.5.
	int minSize = 32;
	int scales = 3;
	//double octaveSigma = sqrt(pow(2*startScale, 2) + pow(preSigma, 2));
	double octaveSigma = 1.6;

	int border = 10;
	int maxSteps = 4;
	//double contrastThresh = 0.03;
	double edgeRatio = 10.0;

	int binsNum = 36;
	double peakThresh = 0.8;

	int gridDim = 4;
	int dirNum = 8;
	int gridSpace = 4;
	double illuThresh = 0.2;

	/////
	//start = clock();

#ifdef DEBUG
	cout << "preblur image." << endl;
#endif

	// modified by styx. (Dec 12, 2007)
	/////////////////////////////////////////////////////////////
	Image image = GaussConvol(preImage, preSigma);
	///////////////////////////////////////////////////////////////


	/////
	//finish = clock();
	//time = (double)(finish - start) / CLOCKS_PER_SEC;
	//cout << "preblur image uses " << time << " s." << endl;

	/////
	//start = clock();

	// Building scale space.
	OctavePyramid pyr;
	//pyr.BuildOctaves(image, minSize, startScale, scales, octaveSigma);
	pyr.BuildOctaves(image, minSize, startScale, scales, octaveSigma);

#ifdef DEBUG
	cout << "There are " << pyr.GetCount() << " octaves." << endl;
#endif

	/////
	//finish = clock();
	//time = (double)(finish - start) / CLOCKS_PER_SEC;
	//cout << "building scale space uses " << time << " s." << endl;

	/////
	//start = clock();

	// Extracting SIFT keypoints in the scale space.
	vector<Keypoint> keys = 
		pyr.BuildKeyList(scales, octaveSigma, 
		border, maxSteps, contrastThresh, edgeRatio, 
		binsNum, peakThresh, 
		gridDim, dirNum, gridSpace, illuThresh);
	cout << "found " << keys.size() << " keypoints in all." << endl;

	/////
	//finish = clock();
	//time = (double)(finish - start) / CLOCKS_PER_SEC;
	//cout << "extracting SIFT keypoints uses " << time << " s." << endl;

	return keys;
}
