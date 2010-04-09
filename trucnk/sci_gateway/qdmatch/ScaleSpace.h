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

// ScaleSpace.h -- Declaration file.
// These classes include detecting SIFT keypoints in the scale-space, 
// assigning orientations and computing descriptor for each keypoint.

#ifndef SCALESPACE_H
#define SCALESPACE_H


#include "Image.h"
#include "GaussConvol.h"
#include "Keypoint.h"
#include "SiftUtil.h"


const double INFMIN = 0.0000001;


/*-------------------------------Each octave---------------------------------*/

class DScaleSpace
{
public:
	DScaleSpace(void);
	~DScaleSpace(void);

	void BuildGaussImgs(Image &first, double firstScale,
		int scales, double sigma);
	void BuildDiffImgs(void);

	DScaleSpace* GetDown(void);
	DScaleSpace* GetUp(void);
	void SetDown(DScaleSpace *dsp);
	void SetUp(DScaleSpace *dsp);

	Image & GetGaussImg(int i);

	Image& operator[](int i);
	const Image& operator[] (int i) const;

	Image & GetLastGaussImg(void);
	int GetCount(void) const;

	vector<Keypoint> FindExtrema(int border);
	vector<Keypoint> FilterExtrema(vector<Keypoint> &keys, 
		int maxSteps, double contrastThresh, double edgeRatio);

	vector<Keypoint> GenerateKeys(int scales, double octaveSigma, 
		int border, int maxSteps, double contrastThresh, double edgeRatio);

	void GenerateDescs(vector<Keypoint> &keys, 
		const vector<Keypoint> &posKeys, 
		int scales, double octaveSigma, int binsNum, double peakThresh, 
		int gridDim, int dirNum, int gridSpace, double illuThresh);

private:
	DScaleSpace(const DScaleSpace &other);
	DScaleSpace & operator =(const DScaleSpace &other);

	void FindLevelExtrema(vector<Keypoint> &keys, 
		Image &below, Image &above, Image &current, 
		int level, int border);
	//bool CheckMin(Image &layer, double c, int x, int y);
	//bool CheckMax(Image &layer, double c, int x, int y);
	void CheckMinMax(Image &layer, double c, int x, int y, 
		int &isMin, int &isMax, bool cLayer);
	bool IsEdge(const Keypoint &key, double edgeRatio);
	bool IsLowContrast(Keypoint &key, double contrastThresh, 
		int steps, int **processed);
	void CalAdjustment(double &xAdjust, double &yAdjust, double &lAdjust, 
		double &dAdjust, const Keypoint &key);
	void SolveLinear(double b[3], double H[3][3], int dim);

	void CalGradImgs(void);
	void DelGradImgs(void);

	vector<Keypoint> AssignOrien(const Keypoint &key, 
		int binsNum, double peakThresh, int scales, double octaveSigma);
	//void AverageWeakBins(double *bins, int binsNum);
	bool ParabolaInter(double &peakPos, double &peakVal, 
		double left, double middle, double right);

	void CreateDescriptor(Keypoint &key, int gridDim, int dirNum, 
		int gridSpace, double illuThresh, int scales, double octaveSigma);
	double BlinearInter(double x, double y, Image &img);
	void ThreshNorm(double *desc, double illuThresh, int dim);

private:
	DScaleSpace *m_down; // pointer to previous octave in octaves list.
	DScaleSpace *m_up; // pointer to the next octave in octaves list.

	Image *m_gaussOctave; // Gaussian images of each octave.
	Image *m_diffOctave; // DoG images of each octave.
	Image m_lastGaussImg; // The last image in Gaussian octave used to build 
	                      // the first image in the next octave.

	double m_baseScale; // subsampling scale for downscaled space.

	int m_count; // length of each DoG octave.

	Image *m_gradMagni; // magnitudes of gradients.
	Image *m_gradOrien; // orientations of gradients.
};


/*--------------------------------Scale space--------------------------------*/

class OctavePyramid
{
public:
	OctavePyramid(void);
	~OctavePyramid(void);

	void BuildOctaves(Image &srcImg, int minSize, double scale, 
		int levelsPerOctave, double octaveSigma);
	DScaleSpace& operator [](int i);
	int GetCount(void) const;

	vector<Keypoint> BuildKeyList(int scales, double octaveSigma, 
		int border, int maxSteps, double contrastThresh, double edgeRatio, 
		int binsNum, double peakThresh, 
		int gridDim, int dirNum, int gridSpace, double illuThresh);

private:
	OctavePyramid(const OctavePyramid &other);
	OctavePyramid &operator =(const OctavePyramid& other);

private:
	DScaleSpace *m_octaves; // pointer to octaves arraylist.
	int m_count; // number of built octaves.
};


#endif
