//////////////////////////////////////////////////////////////////////////////////////////////////
// Author:		Styx(Zhenhui Xu)
// Version:		0.1
// Date:		May 20, 2009 
// Description: Scalespace Defination file.
// Others:      These classes include detecting SIFT keypoints in the scale-space, 
//              assigning orientations and computing descriptor for each keypoint.  
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


#include "ScaleSpace.h"
//#define DEBUG

/*----------------------------------Each octave------------------------------*/

DScaleSpace::DScaleSpace(void)
{
	m_down = 0;
	m_up = 0;

	m_gaussOctave = 0;
	m_diffOctave = 0;
}

DScaleSpace::~DScaleSpace(void)
{
	if (m_gaussOctave != 0)
	{
		delete []m_gaussOctave;
	}
	if (m_diffOctave != 0)
	{
		delete []m_diffOctave;
	}
}

/********************************************************************
   'first' -- the first image over which this octave is built.
   'firstScale' -- downscaled scale.
   'scales' -- scale interval = pow(2.0, 1.0 / scales).
   'sigma' -- initial Gaussian blur sigma.
   This function computes Gaussian Images for this octave.
   We need one more Gaussian image than the number of DoG imags.
   In addition, we need two more for extrema pixel search.
   So we will build ('scales' + 3) Gaussian images by incrementally
   blur the input 'first' image.
********************************************************************/
void DScaleSpace::BuildGaussImgs(Image &first, double firstScale,
								 int scales, double sigma)
{
	// Image first(GaussConvol(srcImg, sigma));

	m_baseScale = firstScale;

	int length = scales + 3;
	m_gaussOctave = new Image[length];
	if (m_gaussOctave == NULL)
	{
		FatalError("BuildGaussImgs -- Allocating memory fails!");
	}

#ifdef DEBUG
	cout << "building level 0 : sigma = " 
		<< (sigma * m_baseScale) << " ..." << endl;
#endif

	// Image prev(first);
	m_gaussOctave[0] = first;

	double k = pow(2, 1.0 / scales);
	double factor = sqrt(k * k - 1);
	double kSigma = sigma;

	for (int i=1; i<length; ++i)
	{
		double gSigma = factor * kSigma;
		kSigma *= k;

#ifdef DEBUG
		cout << "building level " << i << " : sigma = " << 
			(kSigma * m_baseScale) << " ..." << endl;
#endif

		/////
		//clock_t start= clock();

		//prev = GaussConvol(prev, gSigma);
		//m_gaussOctave[i] = prev;
		
		m_gaussOctave[i] = GaussConvol(m_gaussOctave[i-1], gSigma);
		//m_gaussOctave[i] = GaussConvol(GetGaussImg(i-1), gSigma);

		 // m_gaussOctave[i] = GaussConvol(srcImg, kSigma);

		/////
		//clock_t finish = clock();
		//double time = (double)(finish - start) / CLOCKS_PER_SEC;
		//cout << "gauss uses " << time << " s." << endl;
	}

	m_lastGaussImg = m_gaussOctave[length - 3];

	m_count = length - 1;
}

/********************************************************************
   Compute DoG (Difference-of-Gaussian) images for this octave 
   by subtracting two nearby Gaussian images.
********************************************************************/
void DScaleSpace::BuildDiffImgs(void)
{
	m_diffOctave = new Image[m_count];
	if (m_diffOctave == NULL)
	{
		FatalError("BuildDiffImgs -- Allocating memory fails!");
	}

	for (int i=0; i<m_count; ++i)
	{
		m_diffOctave[i] = m_gaussOctave[i+1] - m_gaussOctave[i];
	}
}

DScaleSpace * DScaleSpace::GetDown(void)
{
	return m_down;
}

DScaleSpace * DScaleSpace::GetUp(void)
{
	return m_up;
}

void DScaleSpace::SetDown(DScaleSpace *dsp)
{
	m_down = dsp;
}

void DScaleSpace::SetUp(DScaleSpace *dsp)
{
	m_up = dsp;
}

Image & DScaleSpace::GetGaussImg(int i)
{
	return m_gaussOctave[i];
}

Image& DScaleSpace::operator[](int i)
{
	return m_diffOctave[i];
}

const Image& DScaleSpace::operator[] (int i) const
{
	return m_diffOctave[i];
}

Image & DScaleSpace::GetLastGaussImg(void)
{
	return m_lastGaussImg;
}

int DScaleSpace::GetCount(void) const
{
	return m_count;
}

/********************************************************************
   Detect extrema from DoG images of this octave by comparing every
   pixel with its 26 neighbors, including 8 neighbors in the current 
   image and 9 neighbors in the scale above and below.
********************************************************************/
vector<Keypoint> DScaleSpace::FindExtrema(int border)
{
	vector<Keypoint> keys;
	////////////////////////////////////////////////////////
	Image& tmp = m_diffOctave[0];
	keys.reserve(tmp.GetXDim()*tmp.GetYDim()/20);
	// keys.reserve(1000);
	////////////////////////////////////////////////////////
	for (int level=1; level<(m_count-1); ++level)
	{
		this->FindLevelExtrema(keys, 
			m_diffOctave[level - 1], m_diffOctave[level + 1], 
			m_diffOctave[level], level, border);
	}

	return keys;
}

/********************************************************************
   Find extrema in 'current' image with 'below' and 'above'.
   Ignore extrema out 'border' boundary of the image.
********************************************************************/
void DScaleSpace::FindLevelExtrema(vector<Keypoint> &keys, 
								   Image &below, Image &above, 
								   Image &current, int level, int border)
{
	int xDim = current.GetXDim();
	int yDim = current.GetYDim();

	//assert(border > 0 && border < (Min(xDim, yDim)/2));
	assert(border > 0 && border < (Minimun(xDim, yDim)/2));
	
	for (int y=border; y<(yDim-border); ++y)
	{
		for (int x=border; x<(xDim-border); ++x)
		{
			double c = current(x, y);

			int isMin = 1; // 1 means c is minimum, 0 means not.
			int isMax = 1; // 1 means c is maximum, 0 means not.
		
			this->CheckMinMax(current, c, x, y, isMin, isMax, true);
			this->CheckMinMax(below, c, x, y, isMin, isMax, false);
			this->CheckMinMax(above, c, x, y, isMin, isMax, false);

			if (isMin == 0 && isMax == 0)
			{
				continue;
			}

			/*
			bool isMin = this->CheckMin(current, c, x, y);
			bool isMax = this->CheckMax(current, c, x, y);

			if (isMin)
			{
				isMin = CheckMin(below, c, x, y);
				if (isMin)
				{
					isMin = CheckMin(above, c, x, y);
				}
			}
			else if (isMax)
			{	
				isMax = CheckMax(below, c, x, y);
				if (isMax)
				{
					isMax = CheckMax(above, c, x, y);
				}
			}
			else
			{
				continue;
			}

			if (!(isMin || isMax))
			{
				continue;
			}
			*/

			Keypoint key(x, y, m_baseScale, level);
			keys.push_back(key);
		}
	}
}

/********************************************************************
   Check if a pixel ('x', 'y') with value 'c' is minimum or maximum 
   in the 'layer' image. 
   If 'cLayer' is true, the 'layer' image is the center image, 
   otherwise, the 'layer' image is the above or below image.
********************************************************************/
/*
bool DScaleSpace::CheckMin(Image &layer, double c, int x, int y)
{
	for (int i=-1; i<=1; ++i)
	{
		for (int j=-1; j<=1; ++j)
		{
			if (layer(x + i, y + j) < c)
			{
				return false;
			}
		}
	}

	return true;
}

bool DScaleSpace::CheckMax(Image &layer, double c, int x, int y)
{
	for (int i=-1; i<=1; ++i)
	{
		for (int j=-1; j<=1; ++j)
		{
			if (layer(x + i, y + j) > c)
			{
				return false;
			}
		}
	}

	return true;
}
*/
void DScaleSpace::CheckMinMax(Image &layer, double c, int x, int y, 
							  int &isMin, int &isMax, bool cLayer)
{
	if (isMin == 1)
	{
		if (layer(x - 1, y - 1) <= c ||
			layer(x, y - 1) <= c ||
			layer(x + 1, y - 1) <= c ||
			layer(x - 1, y) <= c ||
			(cLayer ? false : (layer(x, y) < c)) ||
			layer(x + 1, y) <= c ||
			layer(x - 1, y + 1) <= c ||
			layer(x, y + 1) <= c ||
			layer(x + 1, y + 1) <= c)
		{
			isMin = 0;
		}
	}

	if (isMax == 1)
	{
		if (layer(x - 1, y - 1) >= c ||
			layer(x, y - 1) >= c ||
			layer(x + 1, y - 1) >= c ||
			layer(x - 1, y) >= c ||
			(cLayer ? false : (layer(x, y) > c)) ||
			layer(x + 1, y) >= c ||
			layer(x - 1, y + 1) >= c ||
			layer(x, y + 1) >= c ||
			layer(x + 1, y + 1) >= c)
		{
			isMax = 0;
		}
	}
}

/********************************************************************
   'maxSteps' -- the largest times of adjustment loop.
   'contrastThresh' -- If contrast is lower than this, 
                       the keypoint will be discarded.
   'edgeRatio' -- the threshold to discard keypoints along the edge.
   This function computes accurate candidate keypoint localization,  
   then discards unstable keypoints.
********************************************************************/
vector<Keypoint> 
DScaleSpace::FilterExtrema(vector<Keypoint> &keys, int maxSteps, 
						   double contrastThresh, double edgeRatio)
{
	int xDim = m_diffOctave[0].GetXDim();
	int yDim = m_diffOctave[0].GetYDim();
	int x, y;
	// flags to point out duplicates in the same position of one octave.
	int **processed = new int *[yDim];
	if (processed == NULL)
	{
		FatalError("FilterExtrema -- Allocating memory fails!");
	}
	for (y=0; y < yDim; ++y)
	{
		processed[y] = new int[xDim];
		if (processed[y] == NULL)
		{
			FatalError("FilterExtrema -- Allocating memory fails!");
		}
	}

	for (y=0; y<yDim; ++y)
	{
		for(x=0; x<xDim; ++x)
		{
			processed[y][x] = 0;
		}
	}

	// Search all keypoints.
	// Accurate keypoint localization and discard unstable keypoints.
	vector<Keypoint> newKeys;
	////////////////////////////////////////////////////
	newKeys.reserve(keys.size()/3);
	////////////////////////////////////////////////////////

	int num = keys.size();

	for (int i=0; i<num; ++i)
	{
		Keypoint &key = keys[i];

		// Discard keypoints with low contrast.
		if (this->IsLowContrast(key, contrastThresh, maxSteps, processed))
		{
			continue;
		}

		// Discard keypoints along the edge.
		if (this->IsEdge(key, edgeRatio))
		{
			continue;
		}

		newKeys.push_back(key);
	}

	for (y=0; y<yDim; ++y)
	{
		delete [] processed[y];
	}
	delete []processed;

	return newKeys;
}

/********************************************************************
   Determine whether the keypoint is along the edge.
********************************************************************/
bool DScaleSpace::IsEdge(const Keypoint &key, double edgeRatio)
{
	int x = (int)(key.m_x + 0.5);
	int y = (int)(key.m_y + 0.5);
	int level = (int)(key.m_level + 0.5);

	/*
	if (level < 0 || level >= m_count)
	{
		return true;
	}
	if ((x <= 0 || x >= (m_diffOctave[level].GetXDim() - 1)) || 
		(y <= 0 || y >= (m_diffOctave[level].GetYDim() - 1)))
	{
		return true;
	}
	*/
	assert(level > 0 && level < (m_count - 1)); 
	assert(x > 0 && x < (m_diffOctave[level].GetXDim() - 1) && 
		   y > 0 && y < (m_diffOctave[level].GetYDim() - 1));

	Image &current = m_diffOctave[level];

	double dxx = current(x + 1, y) - 2.0 * current(x, y) + current(x - 1, y);
	double dyy = current(x, y + 1) - 2.0 * current(x, y) + current(x, y - 1);
	double dxy = 0.25 * ((current(x + 1, y + 1) - current(x + 1, y - 1)) - 
		(current(x - 1, y + 1) - current(x - 1, y - 1)));

	double trHSq = (dxx + dyy) * (dxx + dyy);
	double detH = dxx * dyy - dxy * dxy;

	if (detH <= 0.0)
	{
		return true;
	}
	else
	{
		double ratio = (edgeRatio + 1.0) * (edgeRatio + 1.0) / edgeRatio;
		if ((trHSq / detH) >= ratio)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
}

/********************************************************************
   Accurate the keypoint localization by looping local adjustment, 
   then determine whether the keypoint is low contrast.
********************************************************************/
bool DScaleSpace::IsLowContrast(Keypoint &key, double contrastThresh, 
								int steps, int **processed)
{
	double dAdjust = 0.0; // adjustment for D(x) value.

	int adjusted = steps; // the largest times of adjustment loop.
	bool needAdjust = true; // If adjusted offset is larger than 0.5, 
	                        // need to adjust again.
	while (needAdjust)
	{
		int x = (int)key.m_x;
		int y = (int)key.m_y;
		int level = (int)key.m_level;

		if (level <= 0 || level >= (m_count - 1))
		{
			return true;
		}
		if ((x <= 0 || x >= (m_diffOctave[level].GetXDim() - 1)) || 
			(y <= 0 || y >= (m_diffOctave[level].GetYDim() - 1)))
		{
			return true;
		}

		double xAdjust = 0.0; // adjustment for x.
		double yAdjust = 0.0; // adjustment for y.
		double lAdjust = 0.0; // adjustment for level.

		this->CalAdjustment(xAdjust, yAdjust, lAdjust, dAdjust, key);

		// If adjusted offset is larger than 0.5, need to adjust again.
		if (fabs(xAdjust) > 0.5 || fabs(yAdjust) > 0.5)
			// || fabs(lAdjust) > 0.5)
		{
			// Already adjusted enough times, give up.
			if (adjusted == 0)
			{
				return true;
			}

			--adjusted;

			// Check that adjustment is less than one pixel step, 
			// otherwise discard the point.
			double distSq = xAdjust * xAdjust + yAdjust * yAdjust;
			if (distSq > 2.0)
			{
				return true;
			}

			// The extremum lies closer to a different sample point, 
			// so change the sample point.
			key.m_x = (int)(key.m_x + xAdjust + 0.5);
			key.m_y = (int)(key.m_y + yAdjust + 0.5);
			//key.SetLevel((int)(key.GetLevel() + lAdjust + 0.5));

			continue;
		}

		if (fabs(lAdjust) >= 0.5)
		{
			return true;
		}
		// Adjusting is over, add adjustment to original value.
		key.m_x = key.m_x + xAdjust;
		key.m_y = key.m_y + yAdjust;
		key.m_level = key.m_level + lAdjust;

		needAdjust = false;
	}

	int x = (int)(key.m_x + 0.5);
	int y = (int)(key.m_y + 0.5);
	int level = (int)(key.m_level + 0.5);

	// We already have a keypoint in this octave for this pixel position.
	if (processed[y][x] != 0)
	{
		return true;
	}

	processed[y][x] = 1;

	// Calculate adjusted D(x) value.
	double d = m_diffOctave[level](x, y) + 0.5 * dAdjust;

	if (fabs(d) < contrastThresh) // Contrast is low.
	{
		return true;
	}
	return false;
}

/********************************************************************
   Calculate adjustment (scaleLevel, y, x) and D(x), respectively 
   return to 'lAdjust', 'yAdjust', 'xAdjust' and 'dAdjust'.
********************************************************************/
void DScaleSpace::CalAdjustment(double &xAdjust, double &yAdjust, 
								double &lAdjust, double &dAdjust, 
								const Keypoint &key)
{
	int x = (int)key.m_x;
	int y = (int)key.m_y;
	int i, j;
	int level = (int)key.m_level;

	assert(level > 0 && level < (m_count - 1) && 
		x > 0 && x < (m_diffOctave[level].GetXDim() - 1) && 
		y > 0 && y < (m_diffOctave[level].GetYDim() - 1));

	double H[3][3];

	Image &below = m_diffOctave[level - 1];
	Image &current = m_diffOctave[level];
	Image &above = m_diffOctave[level + 1];

	H[0][0] = above(x, y) - 2.0 * current(x, y) + below(x, y);
	H[0][1] = H[1][0] = 0.25 * ((above(x, y + 1) - below(x, y + 1)) - 
		(above(x, y - 1) - below(x, y - 1)));
	H[0][2] = H[2][0] = 0.25 * ((above(x + 1, y) - below(x + 1, y)) - 
		(above(x - 1, y) - below(x - 1, y)));
	H[1][1] = current(x, y + 1) - 2.0 * current(x, y) + current(x, y - 1);
	H[1][2] = H[2][1] = 0.25 * 
		((current(x + 1, y + 1) - current(x + 1, y - 1)) - 
		 (current(x - 1, y + 1) - current(x - 1, y - 1)));
	H[2][2] = current(x + 1, y) - 2.0 * current(x, y) + current(x - 1, y);

	double d[3];
	d[0] = 0.5 * (above(x, y) - below(x, y));
	d[1] = 0.5 * (current(x, y + 1) - current(x, y - 1));
	d[2] = 0.5 * (current(x + 1, y) - current(x - 1, y));

	double b[3];
	for (int i=0; i<3; ++i)
	{
		b[i] = - d[i];
	}

	this->SolveLinear(b, H, 3);

	xAdjust = b[2];
	yAdjust = b[1];
	lAdjust = b[0];

	dAdjust = 0.0;
	for (i=0; i<3; ++i)
	{
		dAdjust += d[i] * b[i];
	}
}

/********************************************************************
   Solve linear equations 'H' x = 'b' by Gaussian-elimination method, 
   the solution is also returned to 'b'.
********************************************************************/
void DScaleSpace::SolveLinear(double b[3], double H[3][3], int dim)
{
	int i, j;
	for (j=0; j<(dim-1); ++j)
	{
		int maxIndex = j;
		double maxValue = fabs(H[j][j]);

		for (i=j; i<dim; ++i)
		{
			if (fabs(H[i][j]) > maxValue)
			{
				maxIndex = i;
				maxValue = fabs(H[i][j]);
			}
		}

		if (j != maxIndex)
		{
			for (i=0; i<dim; ++i)
			{
				Swap(H[j][i], H[maxIndex][i]);
			}

			Swap(b[j], b[maxIndex]);
		}

		for (i=j+1; i<dim; ++i)
		{
			double factor = H[i][j] / H[j][j];

			for (int k=0; k<dim; ++k)
			{
				H[i][k] -= factor * H[j][k];
			}

			b[i] -= factor * b[j];
		}
	}

	for (j=dim-1; j>=0; --j)
	{
		double sol = b[j];

		for (int i=dim-1; i>j; --i)
		{ 
			sol -= H[j][i] * b[i];
		}

		b[j] = sol / H[j][j];
	}
}

/********************************************************************
   Calculate gradient magnitudes and orientations of Gaussian images.
********************************************************************/
void DScaleSpace::CalGradImgs(void)
{
	m_gradMagni = new Image[m_count + 1];
	m_gradOrien = new Image[m_count + 1];
	
	for (int i=0; i<(m_count+1); ++i)
	{
		/////
		//clock_t start, finish;
		//double time;

		Image &gaussImg = m_gaussOctave[i];
		Image &magniImg = m_gradMagni[i];
		Image &orienImg = m_gradOrien[i];

		int xDim = gaussImg.GetXDim();
		int yDim = gaussImg.GetYDim();

		/////
		//start = clock();

		magniImg.Allocate(xDim, yDim);
		orienImg.Allocate(xDim, yDim);

		/////
		//finish = clock();
		//time = (double)(finish - start) / CLOCKS_PER_SEC;
		//cout << "allocate uses " << time << " s." << endl;

		/////
		//start = clock();

		for (int y=1; y<(yDim-1); ++y)
		{
			for (int x=1; x<(xDim-1); ++x)
			{
				double diffX = 
					gaussImg(x + 1, y) - gaussImg(x - 1, y);
				double diffY = 
					gaussImg(x, y + 1) - gaussImg(x, y - 1);

				magniImg(x, y) = sqrt(diffX * diffX + diffY * diffY);

				orienImg(x, y) = atan2(diffY, diffX);
				if (fabs(orienImg(x, y) - PI) < INFMIN)
				{
					orienImg(x, y) = -PI;
				}
			}
		}

		/////
		//finish = clock();
		//time = (double)(finish - start) / CLOCKS_PER_SEC;
		//cout << "gradient uses " << time << " s." << endl;
	}
}

/********************************************************************
   Delete gradient magnitude and orientation images.
********************************************************************/
void DScaleSpace::DelGradImgs(void)
{
	if (m_gradMagni != NULL)
	{
		delete []m_gradMagni;
	}
	if (m_gradOrien != NULL)
	{
		delete []m_gradOrien;
	}
}

/********************************************************************
   'key' -- the keypoint to be assigned orientations.
   'binsNum' -- bins number of the orientation histogram.
   'peakThresh' -- threshold with which we find other local peaks.
   'scales' -- frequency of sampling  in scale.
   'octaveSigma' -- initial Gaussian blur sigma.
   This function assigns orientations (between [-PI, PI)) for one 
   keypoint and returns assigned orientation number.
********************************************************************/
vector<Keypoint> 
DScaleSpace::AssignOrien(const Keypoint &key, 
						 int binsNum, double peakThresh, 
						 int scales, double octaveSigma)
{
	int level = (int)(key.m_level + 0.5);
	Image &magnitude = m_gradMagni[level];
	Image &direction = m_gradOrien[level];

	// Build orientation histogram.

	double keyScale = pow(2, key.m_level / scales) * octaveSigma;

	double sigma = 1.5 * keyScale;
	double sigma2Sq = 1.0 / (2.0 * sigma * sigma);
	int radius = (int)(3.0 * sigma + 0.5);
	int radiusSq = radius * radius;

	int keyX = (int)(key.m_x + 0.5);
	int keyY = (int)(key.m_y + 0.5);

//	int xMin = Max(keyX - radius, 1);
//	int xMax = Min(keyX + radius, magnitude.GetXDim() - 1);
//	int yMin = Max(keyY - radius, 1);
//	int yMax = Min(keyY + radius, magnitude.GetYDim() - 1);
	int xMin = Maximun(keyX - radius, 1);
	int xMax = Minimun(keyX + radius, magnitude.GetXDim() - 1);
	int yMin = Maximun(keyY - radius, 1);
	int yMax = Minimun(keyY + radius, magnitude.GetYDim() - 1);

	int x, y, i, j;

	double *bins = new double[binsNum];
	for (i=0; i<binsNum; ++i)
	{
		bins[i] = 0.0;
	}

	for (y=yMin; y<yMax; ++y)
	{
		for (x=xMin; x<xMax; ++x)
		{
			int relX = x - keyX;
			int relY = y - keyY;

			if (relX * relX + relY * relY > radiusSq)
			{
				continue;
			}

			double gWeight = exp(- (relX * relX + relY * relY) * sigma2Sq);
			//gWeight = gWeight / (2.0 * PI * sigma * sigma);

			double dir = direction(x, y);
			if (dir < -PI)
			{
				dir += 2.0 * PI;
			}
			if (dir >= PI)
			{
				dir -= 2.0 * PI;
			}

			// Calculate weight for orientation.
			double idxDir = (dir + PI) * binsNum / (2.0 * PI);
			int binIdxL = (int)idxDir;
			int binIdxR = (binIdxL + 1) % binsNum;
			double dirWeightL = 1.0 - (idxDir - binIdxL);
			double dirWeightR = idxDir - binIdxL;

			bins[binIdxL] += magnitude(x, y) * gWeight * dirWeightL;
			bins[binIdxR] += magnitude(x, y) * gWeight * dirWeightR;
		}
	}

	// Average orientation bins.
	//AverageWeakBins(bins, binsNum);

	// Detect the highest peak.
	double maxValue = 0.0;
	int maxBin = 0;

	for (i=0; i<binsNum; ++i)
	{
		if (bins[i] > maxValue)
		{
			maxValue = bins[i];
			maxBin = i;
		}
	}

	// Interpolate the peak position.
	double peakPos = 0.0;
	double peakVal = maxValue;

	/////////////////////////////////////////////////////////////////
	// I don't think this block has any sense, so I comment it out.
	// Styx (Dec 13, 2007)
	 //this->ParabolaInter(peakPos, peakVal, 
		//bins[(maxBin == 0) ? (binsNum - 1) : (maxBin - 1)], 
		//bins[maxBin], bins[(maxBin + 1) % binsNum]);
	/////////////////////////////////////////////////////////////////


	// Find any other local peak that is within 
	// 'peakThresh' of the highest peak.
	bool *isKey = new bool[binsNum];

	// isKey[maxBin] = true;
	for (i=0; i<binsNum; ++i)
	{
		//isKey[i] = false;

		//if (i == maxBin)
		//{
		//	isKey[i] = true;
		//	continue;
		//}
		if (bins[i] < (peakThresh * peakVal))
		{
			isKey[i] = false;
			continue;
		}

		// not local peak.
		if (bins[i] <= bins[(i == 0) ? (binsNum - 1) : (i - 1)] ||
			bins[i] <= bins[(i + 1) % binsNum])
		{
			isKey[i] = false;
			continue;
		}
		
		isKey[i] = true;
	}

	// Assign orientations to keys.
	double binLen = 2 * PI / binsNum;
	vector<Keypoint> keys;

	for (i=0; i<binsNum; ++i)
	{
		if (!isKey[i])
		{
			continue;
		}

		peakPos = 0.0;
		peakVal = bins[i];

		this->ParabolaInter(peakPos, peakVal, 
			bins[(i == 0) ? (binsNum - 1) : (i - 1)], 
			bins[i], bins[(i + 1) % binsNum]);

		assert(peakPos >= -0.5 && peakPos <= 0.5);

		double keyOrien = (i + peakPos) * binLen - PI;

		if (keyOrien < -PI)
		{
			keyOrien += 2.0 * PI;
		}
		else if (keyOrien >= PI)
		{
			keyOrien -= 2.0 * PI;
		}

		Keypoint newKey(key.m_x, key.m_y, key.m_baseScale, key.m_level);
		newKey.m_scale = keyScale; // * key.m_baseScale;
		newKey.m_orien = keyOrien;
		keys.push_back(newKey);
	}

	delete []isKey;
	delete []bins;

	return keys;
}

/********************************************************************
   Average orientation bins three by three neighbors.
********************************************************************/
/*
void DScaleSpace::AverageWeakBins(double *bins, int binsNum)
{
	for (int n=0; n<4; ++n)
	{
		double first = bins[0];
		double last = bins[binsNum - 1];

		for (int i=0; i<binsNum; ++i)
		{
			double cur = bins[i];
			double next = 
				(i == (binsNum - 1)) ? first : bins[(i + 1) % binsNum];

			bins[i] = (last + cur + next) / 3.0;
			last = cur;
		}
	}
}
*/

/********************************************************************
   Fit a parabola to the three points (-1.0; left), (0.0; middle) 
   and (1.0; right).
   Formula : f(x) = a (x - c)^2 + b.
   where c is the peak offset, b is the peak value.
   If the parabola interpolating successes, return true, 
   otherwise return false.
********************************************************************/
bool DScaleSpace::ParabolaInter(double &peakPos, double &peakVal, 
								double left, double middle, double right)
{
	double a = ((left + right) - 2.0 * middle) / 2.0;

	// not a parabola, a horizontal line.
	if (a == 0.0)
	{
		return false;
	}

	double c = (((left - middle) / a) - 1.0) / 2.0;
	double b = middle - c * c * a;

	// 'middle' is not a peak.
	if (c < -0.5 || c > 0.5)
	{
		return false;
	}

	peakPos = c;
	peakVal = b;

	return true;
}

/********************************************************************
   'key' -- the keypoint for which we create descriptor.
   'gridDim' -- grid dimension of the descriptor, 
                the recommended value is 4.
   'dirNum' -- count of discretized direction, 
               the recommended value is 8.
   'gridSpace' -- grid spacing of the descriptor, 
                  the recommended value is 4.
   'illuThresh' -- threshold to avoid illumination change, 
                   the recommended value is 0.2.
   'scales' -- frequency of sampling  in scale.
   'octaveSigma' -- initial Gaussian blur sigma.
   This function create local image descriptor for the keypoint.
   We build orientation histograms with the gradient samples around 
   the keypoint, then create descriptor with these histogram bins value.
********************************************************************/
void DScaleSpace::CreateDescriptor(Keypoint &key, int gridDim, int dirNum, 
								   int gridSpace, double illuThresh, 
								   int scales, double octaveSigma)
{
	int dim = gridDim * gridDim * dirNum;
	double *desc = new double[dim];
	int i;

	for (i=0; i<dim; ++i)
	{
		desc[i] = 0.0;
	}

	int level = (int)(key.m_level + 0.5);
	Image &magnitude = m_gradMagni[level];
	Image &direction = m_gradOrien[level];
	int xDim = magnitude.GetXDim();
	int yDim = magnitude.GetYDim();

	double xKey = key.m_x;
	double yKey = key.m_y;
	double angle = key.m_orien;

	double dirSpace = 2.0 * PI / dirNum;
	int descWindow = gridDim * gridSpace;
	int radius = descWindow / 2;
	double scaleFactor = key.m_scale / octaveSigma;

	// Gaussian weight sigma and radius.
	double gSigma = descWindow / 2;

	// Search all sample points around the keypoint to create descriptor.
	for (int x=-radius; x<radius; ++x)
	{
		for (int y=-radius; y<radius; ++y)
		{
			// The keypoint is center (0, 0).
			double xS = x + 0.5;
			double yS = y + 0.5;

			double gWeight = exp(- (xS * xS + yS * yS) / 
				(2.0 * gSigma * gSigma));
			//gWeight /= (2.0 * PI * gSigma * gSigma);

			// The coordinates are rotated by 'angle'.
			double xSR = xS * cos(angle) - yS * sin(angle);
			double ySR = xS * sin(angle) + yS * cos(angle);

			// Interpolate magnitude and direction pixel value using bilinear.
			double curX = xKey + xSR * scaleFactor;
			double curY = yKey + ySR * scaleFactor;
			if (curX <= 1 || curX >= (xDim - 2) ||
				curY <= 1 || curY >= (yDim - 2))
			{
				continue;
			}
			double mag = this->BlinearInter(curX, curY, magnitude);
			double ori = this->BlinearInter(curX, curY, direction);
			
			double magW = mag * gWeight;

			// We distribute the value of each gradient sample into 
			// adjacent 8 histogram bins.
			int xIdx[2];
			int yIdx[2];
			int dirIdx[2];
			double xWeight[2];
			double yWeight[2];
			double dirWeight[2];
			for (i=0; i<2; ++i)
			{
				xIdx[i] = 0;
				yIdx[i] = 0;
				dirIdx[i] = 0;
				xWeight[i] = 0.0;
				yWeight[i] = 0.0;
				dirWeight[i] = 0.0;
			}

			// Calculate weights for x, that is (1.0 - d).
			double idxX = (xS + radius - (gridSpace / 2.0)) / gridSpace;
			if (idxX >= 0)
			{
				xIdx[0] = (int)idxX;
				xWeight[0] = 1.0 - (idxX - xIdx[0]);
			}
			if (idxX < (gridDim - 1))
			{
				xIdx[1] = (int)(idxX + 1.0);
				xWeight[1] = 1.0 - (xIdx[1] - idxX);
			}

			// Calculate weights for y, that is (1.0 - d).
			double idxY = (yS + radius - (gridSpace / 2.0)) / gridSpace;
			if (idxY >= 0)
			{
				yIdx[0] = (int)idxY;
				yWeight[0] = 1.0 - (idxY - yIdx[0]);
			}
			if (idxY < (gridDim - 1))
			{
				yIdx[1] = (int)(idxY + 1.0);
				yWeight[1] = 1.0 - (yIdx[1] - idxY);
			}

			// The direction is rotated by 'angle'.	
			double dir = ori - angle;
			if (dir < -PI)
			{
				dir += 2.0 * PI;
			}
			if (dir >= PI)
			{
				dir -= 2.0 * PI;
			}

			// Calculate weight for orientation, that is (1.0 - d).
			double idxDir = (dir + PI) * dirNum / (2.0 * PI);
			if ((int)idxDir == dirNum)
			{
				idxDir -= dirNum;
			}
			dirIdx[0] = (int)idxDir;
			dirIdx[1] = (dirIdx[0] + 1) % dirNum;
			dirWeight[0] = 1.0 - (idxDir - dirIdx[0]);
			dirWeight[1] = idxDir - dirIdx[0];

			// Build orientation histogram, and create descriptor.
			for (int iy = 0 ; iy < 2 ; ++iy)
			{
				for (int ix = 0 ; ix < 2 ; ++ix)
				{
					for (int id = 0 ; id < 2 ; ++id)
					{
						int idx = (xIdx[ix] * gridDim * dirNum) + 
							(yIdx[iy] * dirNum) + dirIdx[id];
						assert(idx >= 0 && idx < 128);

						desc[idx] += 
							magW * xWeight[ix] * yWeight[iy] * dirWeight[id];
					}
				}
			} // end of for
		}
	}

	// Avoid illumination change.
	ThreshNorm(desc, illuThresh, dim);

	// Convert float descriptor values to uchar format.
	uchar *descUchar = new uchar[dim];
	for (i=0; i<dim; ++i)
	{
		int val = (int)(desc[i] * 255.0 + 0.5);
		assert(val >= 0 && val <= 255);

		descUchar[i] = (uchar)val;
	}

	key.CreateDesc(dim, descUchar);

	delete []desc;
	delete []descUchar;
}

/********************************************************************
   Calculate bilinear interpolation pixel value in the image.
********************************************************************/
double DScaleSpace::BlinearInter(double x, double y, Image &img)
{
	int x1 = (int)x;
	int y1 = (int)y;
	int x2 = x1 + 1;
	int y2 = y1 + 1;

	assert(x1 >= 0 && x2 >= 0);
	assert(x2 < img.GetXDim() && y2 < img.GetYDim());

	double val = 
		(x2 - x) * (y2 - y) * img(x1, y1) + 
		(x - x1) * (y2 - y) * img(x2, y1) + 
		(x2 - x) * (y - y1) * img(x1, y2) + 
		(x - x1) * (y - y1) * img(x2, y2);

	return val;
}

/********************************************************************
   To avoid linear illumination change, we normalize the descriptor.
   To avoid non-linear illumination change, we threshold the value 
   of each descriptor element to 'illuThresh', then normalize again.
********************************************************************/
void DScaleSpace::ThreshNorm(double *desc, double illuThresh, int dim)
{
	// Normalize the descriptor, and threshold 
	// value of each element to 'illuThresh'.

	double norm = 0.0;
	int i;

	for (i=0; i<dim; ++i)
	{
		norm += desc[i] * desc[i];
	}

	norm = sqrt(norm);
	assert(norm != 0);

	for (i=0; i<dim; ++i)
	{
		desc[i] /= norm;

		if (desc[i] > illuThresh)
		{
			desc[i] = illuThresh;
		}
	}

	// Normalize again.

	norm = 0.0;

	for (i=0; i<dim; ++i)
	{
		norm += desc[i] * desc[i];
	}

	norm = sqrt(norm);
	assert(norm != 0);

	for (i=0; i<dim; ++i)
	{
		desc[i] /= norm;
	}
}

/********************************************************************
   Generate keypoints for this octave.
********************************************************************/
vector<Keypoint> 
DScaleSpace::GenerateKeys(int scales, double octaveSigma, 
						  int border, int maxSteps, 
						  double contrastThresh, double edgeRatio)
{
	/////
	//clock_t start, finish;
	//double time;

	/////
	//start = clock();

	// Detect extrema for candidate keypoints.
	vector<Keypoint> initialKeys = this->FindExtrema(border);

	/////
	//finish = clock();
	//time = (double)(finish - start) / CLOCKS_PER_SEC;
	//cout << "finding extrema " << initialKeys.size() <<
	//	" uses " << time << " s." << endl;

	/////
	//start = clock();

	// Accurate keypoint localization and discard unstable keypoints.
	vector<Keypoint> keys = this->FilterExtrema(initialKeys, 
		maxSteps, contrastThresh, edgeRatio);

	/////
	//finish = clock();
	//time = (double)(finish - start) / CLOCKS_PER_SEC;
	//cout << "filtering extrema " << keys.size() << 
	//	" uses " << time << " s." << endl;

	return keys;
}

/********************************************************************
   Generate keypoint descriptors for this octave.
********************************************************************/
void DScaleSpace::GenerateDescs(vector<Keypoint> &keys, 
								const vector<Keypoint> &posKeys, 
								int scales, double octaveSigma, 
								int binsNum, double peakThresh, 
								int gridDim, int dirNum, int gridSpace, 
								double illuThresh)
{
	/////
	//clock_t start, finish;
	//double time;

	/////
	//start = clock();

	// Pre-compute gradient magnitudes and orientations.
	this->CalGradImgs();

	/////
	//finish = clock();
	//time = (double)(finish - start) / CLOCKS_PER_SEC;
	//cout << "computing gradients uses " << time << " s." << endl;

	/////
	//start = clock();

	// Compute descriptors for every keypoint.
	int num = posKeys.size();
	for (int i=0; i<num; ++i)
	{
		// Assign orientations for each keypoint.
		vector<Keypoint> orienKeys = this->AssignOrien(posKeys[i], 
			binsNum, peakThresh, scales, octaveSigma);

		int orienNum = orienKeys.size();
		for (int j=0; j<orienNum; ++j)
		{
			Keypoint &key = orienKeys[j];

			// Create descriptor for each keypoint.
			this->CreateDescriptor(key, 
				gridDim, dirNum, gridSpace, illuThresh, 
				scales, octaveSigma);

			// Adjust position of each keypoint to original 
			// image size by base scale value.
			key.m_x *= key.m_baseScale;
			key.m_y *= key.m_baseScale;
			key.m_scale *= key.m_baseScale;
			
			keys.push_back(key);
		}
	}

	/////
	//finish = clock();
	//time = (double)(finish - start) / CLOCKS_PER_SEC;
	//cout << "creating descriptors uses " << time << " s." << endl;

	// Delete gradient magnitudes and orientations.
	this->DelGradImgs();
}


/*--------------------------------Scale space--------------------------------*/

OctavePyramid::OctavePyramid(void)
{
	m_octaves = 0;
	m_count = 0;
}

OctavePyramid::~OctavePyramid(void)
{
	DScaleSpace *pcur, *pdown;
	pcur = m_octaves;

	while (pcur != 0)
	{
		pdown = pcur->GetDown();

		delete pcur;

		pcur = pdown;
	}
}

/********************************************************************
   'srcImg' -- source image.
   'minSize' -- minimum size of the image in scale space.
   'scale' -- downscaled scale.
   'levelsPerOctave' -- scale interval = pow(2.0, 1/levelsPerOctave).
   'octaveSigma' -- initial Gaussian blur sigma.
   This function builds the largest possible number of octaves.
   Each octave is downscaled by 0.5 and sigma doubled for the next 
   octave until the size of image is less than minSize.
********************************************************************/
void OctavePyramid::BuildOctaves(Image& srcImg, int minSize, double scale, 
								 int levelsPerOctave, double octaveSigma)
{
	m_octaves = NULL;
	m_count = 0;
	Image preImg(srcImg);
	//Image preImg(GaussConvol(srcImg, octaveSigma));

	while (preImg.GetXDim() >= minSize && preImg.GetYDim() >= minSize)
	{
#ifdef DEBUG
		cout << "building octave "<< m_count << "..." << endl;
#endif

		DScaleSpace *dsp = new DScaleSpace;

		if (dsp == 0)
		{
			FatalError("BuildOctaves -- Allocating memory fails!");
		}

		dsp->BuildGaussImgs(preImg, scale, levelsPerOctave, octaveSigma);
		dsp->BuildDiffImgs();

		if (m_octaves != 0)
		{
			m_octaves->SetUp(dsp);
		}
		dsp->SetDown(m_octaves);
		m_octaves = dsp;

		++m_count;

		preImg = dsp->GetLastGaussImg().HalfScale();

		scale *= 2.0;
	}
}

DScaleSpace& OctavePyramid::operator[](int i)
{
	DScaleSpace *octave = m_octaves;
	for (int n=1; n<(m_count-i); ++n)
	{
		octave = octave->GetDown();
	}
	return *octave;
}

int OctavePyramid::GetCount(void) const
{
	return m_count;
}

/********************************************************************
   Build keypoints list by linking all keypoints generated 
   from each octave.
   The parameters are the same as those in class 'DScaleSpace'.
********************************************************************/
vector<Keypoint> 
OctavePyramid::BuildKeyList(int scales, double octaveSigma, 
							int border, int maxSteps, 
							double contrastThresh, double edgeRatio, 
							int binsNum, double peakThresh, 
							int gridDim, int dirNum, int gridSpace, 
							double illuThresh)
{
	/////
	//clock_t start, finish;
	//double time;

	vector<Keypoint> allKeys;
	
	// Modified by Styx, Dec 18, 2007
	// Reserve space for vector (speed up)
	Image& tmp = m_octaves[0].GetGaussImg(0);
	allKeys.reserve(tmp.GetXDim()*tmp.GetYDim()/100);
	//allKeys.reserve(5000);

	for (int i=0; i<m_count; ++i)
	{
		DScaleSpace &octave = (*this)[i];
		
		/////
		//start = clock();

		// generate keypoints.
		vector<Keypoint> octaveKeys = 
			octave.GenerateKeys(scales, octaveSigma, 
			border, maxSteps, contrastThresh, edgeRatio);

		/////
		//finish = clock();
		//time = (double)(finish - start) / CLOCKS_PER_SEC;
		//cout << "generating keypoints uses " << time << " s." << endl;

		/////
		//start = clock();

		// generate descriptors.
#ifdef DEBUG
		int beforeNum = allKeys.size();
#endif
		octave.GenerateDescs(allKeys, octaveKeys, 
			scales, octaveSigma, binsNum, peakThresh, 
			gridDim, dirNum, gridSpace, illuThresh);
#ifdef DEBUG
		int afterNum = allKeys.size();
#endif

		/////
		//finish = clock();
		//time = (double)(finish - start) / CLOCKS_PER_SEC;
		//cout << "generating descriptors uses " << time << " s." << endl;

#ifdef DEBUG
		cout << "found " << (afterNum - beforeNum ) << 
			" keypoints in octave " << i << " ." << endl;
#endif
	}

	return allKeys;
}
