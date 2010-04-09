//////////////////////////////////////////////////////////////////////////////////////////////////
// Author:		Styx(Zhenhui Xu)
// Version:		0.1
// Date:		May 20, 2009 
// Description: This funciton is used to match two set of image points.
// Others:          
// Function List:
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

#include "Correspond.h"
#include "GlobalPara.h"


/********************************************************************
   This function computes correspondence between two image points, 
   using some steps to discard outliers and readjust the distribution.

   matches -- output matches after remove outliers.
   initialMatches -- input initial matches.
   xDim -- width of the first image.
   yDim -- height of the first image.
********************************************************************/
int Correspond(NyxMat<double> &F, vector<Match> &matches, 
			   const vector<Match> &initialMatches, 
			   int xDim, int yDim)
{
	matches.clear();
	int matchNum = 0;

	int i = 0;
	bool isHomog = false;

	// Sort the matches by an increasing order of distances.
	//sort(initialMatches.begin(), initialMatches.end());

	vector<int> inliers;
	inliers.reserve(initialMatches.size());

	NyxMat<double> fund;

	// RANSAC using fundamental matrix model.
	int flag = RansacFitFund(fund, inliers, initialMatches, RANSAC_THRESH);
	if (flag == 0)
	{
		cout << "RANSAC failed." << endl;

		matches = initialMatches;
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				fund(i,j) = 0;
			}
		}
		return matchNum;
	}

	//vector<Match> matches1;
	matchNum = inliers.size();

	// Find more matches. 
	vector<Match> matches2;
	matches2.reserve(initialMatches.size()/2);

	matchNum = EpiLineConstrn(matches2, initialMatches, fund, EPIPOLAR_CONSTRAINT_THRESH);

	// Distribute features across the image using non-minimal suppression.
	vector<Match> matches3;
	matches3.reserve(matchNum);

	int radius = 5;
	matchNum = NonMinSupp(matches3, matches2, xDim, yDim, radius);
	matches2.clear();

	// Remove outliers again.
	vector<Match> matches4;
	matches4.reserve(matchNum);

	double thresh3 = 0.0005;
	flag = RansacFitFund(F, inliers, matches3, thresh3);

	if (flag == 0)
	{
		cout << "RANSAC failed." << endl;
		matches4 = matches3;
	}
	matchNum = EpiLineConstrn(matches4, initialMatches, F, EPIPOLAR_CONSTRAINT_THRESH-0.5);

	matches = matches4;

	return matchNum;
}


/********************************************************************
   Compute the Zero-mean Normalized Cross-Correlation coefficients 
   between the matched window pairs. To get sub-pixel accuracy, 
   we use parabola interpolation to refine the cross-correlation 
   measure values. The final matches are sorted by a decreasing 
   order of correlation coefficients.
********************************************************************/
int CorrMatch(vector<Match> &corrMatches, 
			  const vector<Match> &matches, Image &im1, Image &im2)
{
	corrMatches.clear();

	// NEIGHBOR: the range of neighbors to be searched : (-neighb~neighb)^2.

	// the rectangular correlation window size : window*window.
	int window = ZNCC_WINDOW;
	int radius = (window - 1) / 2;

	// the threshold to remove matches with low correlation scores.
	//double thresh = SCORE_THRESH;

	int xDim1 = im1.GetXDim();
	int yDim1 = im1.GetYDim();
	int xDim2 = im2.GetXDim();
	int yDim2 = im2.GetYDim();

	int matchNum = matches.size();

	// for each match pair.
	for (int i=0; i<matchNum; ++i)
	{
		int x1 = (int)(matches[i].GetTarget().GetX() + 0.5);
		int y1 = (int)(matches[i].GetTarget().GetY() + 0.5);
		int x2 = (int)(matches[i].GetFound().GetX() + 0.5);
		int y2 = (int)(matches[i].GetFound().GetY() + 0.5);

		// Ignore the matched points which cannot get full correlation windows.
		if ((x1-NEIGHBOR+radius) >= xDim1 || (x1+NEIGHBOR-radius) < 0 || 
			(y1-NEIGHBOR+radius) >= yDim1 || (y1+NEIGHBOR-radius) < 0 || 
			(x2-NEIGHBOR+radius) >= xDim2 || (x2+NEIGHBOR-radius) < 0 || 
			(y2-NEIGHBOR+radius) >= yDim2 || (y2+NEIGHBOR-radius) < 0)
		{
			continue;
		}

		// refined coordinate values.
		int adjustX2 = x2;
		int adjustY2 = y2;

		// Store the maximum correlation coefficient.
		double maxCoef = -1.0;

		// for each sample point in the neighborhood.
		for (int x=-NEIGHBOR; x<=NEIGHBOR; ++x)
		{
			for (int y=-NEIGHBOR; y<=NEIGHBOR; ++y)
			{
				// modified by styx, Dec 17, 2007
				// ignore candidate match-points beyond the boundary
				if ( (x2+x) < radius || (x2+x) >= xDim2-radius ||
					(y2+y) < radius || (y2+y) >= yDim2-radius)
					continue;

				double coef = NormCrossCorr(x1, y1, 
					x2+x, y2+y, radius, im1, im2);
				
				// Remove matches with low correlation scores.
				if (coef < SCORE_THRESH)
				{
					continue;
				}
				
				// Store the best match with the largest score.
				if (coef > maxCoef)
				{
					maxCoef = coef;
					adjustX2 = x2 + x;
					adjustY2 = y2 + y;
				}
			}
		}

		if (maxCoef < SCORE_THRESH)
		{
			continue;
		}

		// Refine matches by parabola interpolation.
		double fineX2, fineY2, fineCoef;

		double coef1 = NormCrossCorr(x1, y1, 
			adjustX2-1, adjustY2, radius, im1, im2);
		double coef2 = NormCrossCorr(x1, y1, 
			adjustX2+1, adjustY2, radius, im1, im2);
		if (ParabolaInter(fineX2, fineCoef, coef1, maxCoef, coef2))
		{
			fineX2 += adjustX2;
		}
		else
		{
			//fineX2 = adjustX2;
			continue;
		}
		
		coef1 = NormCrossCorr(x1, y1, 
			adjustX2, adjustY2-1, radius, im1, im2);
		coef2 = NormCrossCorr(x1, y1, 
			adjustX2, adjustY2+1, radius, im1, im2);
		if (ParabolaInter(fineY2, fineCoef, coef1, maxCoef, coef2))
		{
			fineY2 += adjustY2;
		}
		else
		{
			//fineY2 = adjustY2;
			continue;
		}

		Point p1(x1, y1);
		Point p2(fineX2, fineY2);

		// final results.
		// Store correlation coefficients in the matches.
		Match match(p1, p2, maxCoef);
		corrMatches.push_back(match);
	}

	// Sort the matches by a decreasing order of correlation coefficients.
	sort(corrMatches.begin(), corrMatches.end(), less<Match>());
	
	return corrMatches.size();
}

/********************************************************************
   Compute Zero-mean Normalized Cross-Correlation coefficient as : 
                           sum_W( (W1-E(W1))*(W2-E(W2)) ) 
   ZNCC(W1,W2)= ------------------------------------------------------
                sqrt( sum_W1( (W1-E(W1))^2 )* sum_W2( (W2-E(W2))^2 ) )
********************************************************************/
double NormCrossCorr(int x1, int y1, int x2, int y2, int radius, 
					 Image &im1, Image &im2)
{
	int window = 2 * radius + 1;
	int windowSq = window * window;

	int xDim1 = im1.GetXDim();
	int yDim1 = im1.GetYDim();
	int xDim2 = im2.GetXDim();
	int yDim2 = im2.GetYDim();

	//int binNum = windowSq;

	//assert(x1 >= 0 && x1 < xDim1 && y1 >= 0 && y1 < yDim1 && 
	//	   x2 >= 0 && x2 < xDim2 && y2 >= 0 && y2 < yDim2);

	double mean1 = 0.0;
	double mean2 = 0.0;

	if ((x1-radius) < 0 || (x1+radius) >= xDim1 || 
		(y1-radius) < 0 || (y1+radius) >= yDim1 ||
		(x2-radius) < 0 || (x2+radius) >= xDim2 ||
		(y2-radius) < 0 || (y2+radius) >= yDim2)
		return 0;

	for (int x=-radius; x<=radius; ++x)
	{
		for (int y=-radius; y<=radius; ++y)
		{
			//if ((x1+x) < 0 || (x1+x) >= xDim1 || 
			//	(y1+y) < 0 || (y1+y) >= yDim1 ||
			//	(x2+x) < 0 || (x2+x) >= xDim2 ||
			//	(y2+y) < 0 || (y2+y) >= yDim2)
			//{
			//	return 0;
			//	//--binNum;
			//	//continue;
			//}

			mean1 += im1(x1+x, y1+y);
			mean2 += im2(x2+x, y2+y);
		}
	}

	mean1 /= (double)windowSq;
	mean2 /= (double)windowSq;
	//mean1 /= (double)binNum;
	//mean2 /= (double)binNum;

	double diffCross = 0.0;
	double diffSq1 = 0.0;
	double diffSq2 = 0.0;

	double diff1, diff2;

	for (int x=-radius; x<=radius; ++x)
	{
		for (int y=-radius; y<=radius; ++y)
		{
			//if ((x1+x) < 0 || (x1+x) >= xDim1 || 
			//	(y1+y) < 0 || (y1+y) >= yDim1 ||
			//	(x2+x) < 0 || (x2+x) >= xDim2 ||
			//	(y2+y) < 0 || (y2+y) >= yDim2)
			//{
			//	continue;
			//}

			diff1 = im1(x1+x, y1+y) - mean1;
			diff2 = im2(x2+x, y2+y) - mean2;
			
			diffCross += diff1 * diff2;
			diffSq1 += diff1 * diff1;
			diffSq2 += diff2 * diff2;
		}
	}

	double coef = diffCross / sqrt(diffSq1 * diffSq2);

	return fabs(coef);
}

/********************************************************************
   Fit a parabola to the three points (-1.0; left), (0.0; middle) 
   and (1.0; right).
   Formula : f(x) = a (x - c)^2 + b.
   where c is the peak offset, b is the peak value.
   If the parabola interpolating successes, return true, 
   otherwise return false.
********************************************************************/
bool ParabolaInter(double &peakPos, double &peakVal, 
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
   reMatches -- remained matches.
   matches -- original matches.
   xDim -- width of the first image.
   yDim -- height of the first image.
   radius -- radius of region considered in non-minimal suppression.
   This function preforms non-minimal suppression for features.
********************************************************************/
int NonMinSupp(vector<Match> &reMatches, 
			   const vector<Match> &matches, 
			   int xDim, int yDim, int radius)
{
	reMatches.clear();

	// Store the closest distance cost of the matched points.
	Image costImg(xDim, yDim);

	// Initialize.
	for (int y=0; y<yDim; ++y)
	{
		for (int x=0; x<xDim; ++x)
		{
			costImg(x, y) = 255.0;
		}
	}

	int matchNum = matches.size();

	vector<int> ind;
	//////////////////////////////////////////////////////////////////
	ind.reserve(matchNum);
	//////////////////////////////////////////////////////////////////

	for (int i=0; i<matchNum; ++i)
	{
		int x1 = (int)(matches[i].GetTarget().GetX() + 0.5);
		int y1 = (int)(matches[i].GetTarget().GetY() + 0.5);
		//int x2 = (int)(matches[i].GetFound().GetX() + 0.5);
		//int y2 = (int)(matches[i].GetFound().GetY() + 0.5);

		// Exclude points within radius of the image boundary.
		if (x1 >= (xDim-radius) || x1 < radius || 
			y1 >= (yDim-radius) || y1 < radius)
			//x2 >= (xDim-radius) || x2 < radius || 
			//y2 >= (yDim-radius) || y2 < radius)
		{
			continue;
		}

		// Store distance in the image.
		if (costImg(x1, y1)>=(255.0-EPS) && costImg(x1, y1)<=(255.0+EPS))
		{
			costImg(x1, y1) = matches[i].GetDist();

			ind.push_back(i);
		}
	}

	int n = ind.size();

	for (int ii=0; ii<n; ++ii)
	{
		int i = ind[ii];
		
		int x1 = (int)(matches[i].GetTarget().GetX() + 0.5);
		int y1 = (int)(matches[i].GetTarget().GetY() + 0.5);
		//int x2 = (int)(matches[i].GetFound().GetX() + 0.5);
		//int y2 = (int)(matches[i].GetFound().GetY() + 0.5);

		int label = 1;
		
		// non-minimal suppression.
		for (int x=-radius; x<=radius; ++x)
		{
			for (int y=-radius; y<=radius; ++y)
			{
				if (costImg(x1, y1) > costImg(x1+x, y1+y))
				{
					label = 0;
					break;
				}
			}

			if (label == 0)
			{
				break;
			}
		}

		// Remain the minimal match within the radius.
		if (label == 1)
		{
			reMatches.push_back(matches[i]);
		}
	}

	return reMatches.size();
}

/********************************************************************
   reMatches -- remained matches.
   matches -- original matches.
   xDim -- width of the first image.
   yDim -- height of the first image.
   radius -- radius of region considered in non-maximal suppression.
   This function preforms non-maximal suppression for features.
********************************************************************/
int NonMaxSupp(vector<Match> &reMatches, 
			   const vector<Match> &matches, 
			   int xDim, int yDim, int radius)
{
	reMatches.clear();

	// initial cost value.
	double iValue = -1.0;

	// Store the ZNCC correlation cost of the matched points.
	Image costImg(xDim, yDim);

	// Initialize.
	for (int y=0; y<yDim; ++y)
	{
		for (int x=0; x<xDim; ++x)
		{
			costImg(x, y) = iValue;
		}
	}

	int matchNum = matches.size();

	vector<int> ind;
	///////////////////////////////////////////////////////////////
	ind.reserve(matchNum);
	///////////////////////////////////////////////////////////////

	for (int i=0; i<matchNum; ++i)
	{
		int x1 = (int)(matches[i].GetTarget().GetX() + 0.5);
		int y1 = (int)(matches[i].GetTarget().GetY() + 0.5);
		//int x2 = (int)(matches[i].GetFound().GetX() + 0.5);
		//int y2 = (int)(matches[i].GetFound().GetY() + 0.5);

		// Exclude points within radius of the image boundary.
		if (x1 >= (xDim-radius) || x1 < radius || 
			y1 >= (yDim-radius) || y1 < radius)
			//x2 >= (xDim-radius) || x2 < radius || 
			//y2 >= (yDim-radius) || y2 < radius
		{
			continue;
		}

		// Store distance in the image.
		if (costImg(x1, y1) >= (iValue-EPS) && 
			costImg(x1, y1) <= (iValue+EPS))
		{
			costImg(x1, y1) = matches[i].GetDist();

			ind.push_back(i);
		}
	}

	int n = ind.size();

	for (int ii=0; ii<n; ++ii)
	{
		int i = ind[ii];
		
		int x1 = (int)(matches[i].GetTarget().GetX() + 0.5);
		int y1 = (int)(matches[i].GetTarget().GetY() + 0.5);
		//int x2 = (int)(matches[i].GetFound().GetX() + 0.5);
		//int y2 = (int)(matches[i].GetFound().GetY() + 0.5);

		int label = 1;

		// non-maximal suppression.
		for (int x=-radius; x<=radius; ++x)
		{
			for (int y=-radius; y<=radius; ++y)
			{
				if (costImg(x1, y1) < costImg(x1+x, y1+y))
				{
					label = 0;
					break;
				}
			}

			if (label == 0)
			{
				break;
			}
		}

		// Remain the maximal match within the radius.
		if (label == 1)
		{
			reMatches.push_back(matches[i]);
		}
	}

	return reMatches.size();
}



/********************************************************************
   Compute epipolar lines in two images, and find matches with 
   consistent epipolar line constrain.
********************************************************************/
int EpiLineConstrn(vector<Match> &reMatches, 
				   const vector<Match> &matches, 
				   const NyxMat<double> &fund, double thresh)
{
	assert(fund.GetRow() == 3 && fund.GetCol() == 3);

	reMatches.clear();

	int matchNum = matches.size();

	for (int i=0; i<matchNum; ++i)
	{
		if (IsInLimit(matches[i], fund, thresh))
		{
			reMatches.push_back(matches[i]);
		}
	}
	
	return reMatches.size();
}

int EpiLineConstrn(vector<Match> &reMatches, vector<int>& inliers,
				   const vector<Match> &matches, 
				   const NyxMat<double> &fund, double thresh)
{
	assert(fund.GetRow() == 3 && fund.GetCol() == 3);

	reMatches.clear();
	inliers.clear();

	int matchNum = matches.size();

	for (int i=0; i<matchNum; ++i)
	{
		if (IsInLimit(matches[i], fund, thresh))
		{
			reMatches.push_back(matches[i]);
			inliers.push_back(1);
		}
		else
		{
			inliers.push_back(0);
		}
	}
	
	return reMatches.size();
}

// To see whether a match is within the epipolar constaint
bool IsInLimit(const Match& match, const NyxMat<double> &fund, double thresh)
{
	double x1 = match.GetTarget().GetX();
	double y1 = match.GetTarget().GetY();
	double x2 = match.GetFound().GetX();
	double y2 = match.GetFound().GetY();

	// epipolar line in the second image : l2 = F*P1.
	double l2x = fund(0,0) * x1 + fund(0,1) * y1 + fund(0,2);
	double l2y = fund(1,0) * x1 + fund(1,1) * y1 + fund(1,2);
	double l2z = fund(2,0) * x1 + fund(2,1) * y1 + fund(2,2);

	// epipolar line in the first image : l1 = F'*P2.
	double l1x = fund(0,0) * x2 + fund(1,0) * y2 + fund(2,0);
	double l1y = fund(0,1) * x2 + fund(1,1) * y2 + fund(2,1);
	double l1z = fund(0,2) * x2 + fund(1,2) * y2 + fund(2,2);

	// distance error between the point and 
	// the corresponding epipolar line.
	double error1 = (x1*l1x + y1*l1y + l1z) / sqrt(l1x*l1x + l1y*l1y);
	double error2 = (x2*l2x + y2*l2y + l2z) / sqrt(l2x*l2x + l2y*l2y);
	if ((fabs(error1) < thresh) && (fabs(error2) < thresh))
	{
		return true;
	}

	return false;
}

bool IsInLimit(double x1, double y1,
			   double x2, double y2,
			   const NyxMat<double> &fund, double thresh)
{
	// epipolar line in the second image : l2 = F*P1.
	double l2x = fund(0,0) * x1 + fund(0,1) * y1 + fund(0,2);
	double l2y = fund(1,0) * x1 + fund(1,1) * y1 + fund(1,2);
	double l2z = fund(2,0) * x1 + fund(2,1) * y1 + fund(2,2);

	// epipolar line in the first image : l1 = F'*P2.
	double l1x = fund(0,0) * x2 + fund(1,0) * y2 + fund(2,0);
	double l1y = fund(0,1) * x2 + fund(1,1) * y2 + fund(2,1);
	double l1z = fund(0,2) * x2 + fund(1,2) * y2 + fund(2,2);

	// distance error between the point and 
	// the corresponding epipolar line.
	double error1 = (x1*l1x + y1*l1y + l1z) / sqrt(l1x*l1x + l1y*l1y);
	double error2 = (x2*l2x + y2*l2y + l2z) / sqrt(l2x*l2x + l2y*l2y);
	if ((fabs(error1) < thresh) && (fabs(error2) < thresh))
	{
		return true;
	}

	return false;
}
