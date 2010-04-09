//////////////////////////////////////////////////////////////////////////////////////////////////
// Author:		Styx(Zhenhui Xu)
// Version:		0.1
// Date:		May 20, 2009 
// Description: Quasi-dense propagation algorithms.
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


#include "QuasiDense.h"
#include "GlobalPara.h"


/********************************************************************
   fundM -- output, estimated fundamental matrix.
   matches -- output, sorted by a decreasing order of ZNCC score.
   corrMatches -- input matches with ZNCC correlation scores.
   im1 -- the first input image.
   im2 -- the second input image.
   This function computes quasi-dense correspondences between two 
   images given initial matches, it returns the number of matches.
   This algorithm follows the paper by Maxime Lhuillier & Long Quan: 
   "A Quasi-Dense Approach to Surface Reconstruction from 
   Uncalibrated Images" published in TPAMI 2005.
********************************************************************/
int QuasiDense(NyxMat<double> &fundM, vector<Match> &matches, 
			   const vector<Match> &corrMatches, 
			   Image &im1, Image &im2, int times, 
			   vector<Match>* pMatch01, char* match12File)
{
	int i, j;
	fundM = NyxMat<double>(3, 3, 0);

	matches.clear();
	int matchNum = 0;

	// constrained propagation. //////////////////////////////////
	NyxMat<double> fund;
	vector<int> inliers;

	// In order to save time, I first estimate fundamental matrix using 
	// initial matches extracted by SIFT. So, in fact, unconstrained propagation
	// is ignored.
	int flag = RansacFitFund(fund, inliers, corrMatches, F_THRESH);
	if (flag == 0)
	{
		cout << "RANSAC failed." << endl;

		matches = corrMatches;
		return matches.size();
	}

	vector<Match> denseMatches;
	//denseMatches.reserve(Min(im1.GetXDim()*im1.GetYDim(), im2.GetXDim()*im2.GetYDim())/4);
	denseMatches.reserve(Minimun(im1.GetXDim()*im1.GetYDim(), im2.GetXDim()*im2.GetYDim())/4);

	matchNum = Propagate(denseMatches, corrMatches, fund, E_THRESH, im1, im2);


	// Resample quasi-dense pixel correspondences 
	// by local geometric constraints. //////////////////////////////
	vector<Match> sampledMatches;
	sampledMatches.reserve(corrMatches.size()+im1.GetXDim()*im1.GetYDim()/8/8);

	if (pMatch01 == 0)
		matchNum = Resample(sampledMatches, denseMatches, im1, im2, times);
	else
		matchNum = Resample(sampledMatches, denseMatches, 
		im1, im2, times, pMatch01, match12File);		
	
//	printf("Resample done ......\n");

	// Estimate fundamental matrix again. ///////////////////////////
	/////

	denseMatches.clear();
	flag = RansacFitFund(fund, inliers, sampledMatches, F_THRESH);
	if (flag == 0)
	{
		cout << "RANSAC failed." << endl;

		matches = sampledMatches;
		return matches.size();
	}

	//// Display fundamental matrix.
	//cout << "F = " << endl;
	//for (i=0; i<3; ++i)
	//{
	//	for (j=0; j<3; ++j)
	//	{
	//		cout << fund[i][j] << "; ";
	//	}
	//	cout << endl;
	//}

	if (pMatch01 == 0)
	{
		matchNum = EpiLineConstrn(matches, sampledMatches, fund, E_THRESH);
	}
	else
	{
		vector<Match> tmpMatch;
		vector<Match>& refMatch01 = *pMatch01;
		matchNum = EpiLineConstrn(matches, inliers, sampledMatches, fund, E_THRESH);
		for (int i = 0; i < inliers.size(); ++i)
		{
			if (inliers[i] == 1)
			{
				tmpMatch.push_back(refMatch01[i]);
			}
		}
		(*pMatch01) = tmpMatch;
	}
	
	
	fundM = fund;

	//matches = sampledMatches;

	// Sort final matches in a decreasing order of correlation coefficients.
	// sort(matches.begin(), matches.end(), less<Match>());

	return matches.size();
}


/********************************************************************
   Propagate input sparse matches to get quasi-dense matches.
   This algorithm uses best-first strategy, input match seeds 
   must include ZNCC correlation score.
********************************************************************/
int Propagate(vector<Match> &denseMatches, 
			  const vector<Match> &matches, NyxMat<double> &fund, double thresh, 
			  Image &im1, Image &im2)
{
	/////
	//clock_t start, finish;
	//double time;
	
	//int neighb = 2;				// radius of the neighborhood region.
	int neighb = PROPA_NEIGHBOR;
	int radius = PROPA_RADIUS;		// radius of the ZNCC window.
	int dThresh = D_THRESH;			// disparity gradient limit.
	double cThresh = C_THRESH;		// ZNCC value threshold.
	double sThresh = S_THRESH;		// confidence measure threshold.

	int xDim1 = im1.GetXDim();
	int yDim1 = im1.GetYDim();
	int xDim2 = im2.GetXDim();
	int yDim2 = im2.GetYDim();

	// Initialize uniqueness flags for image 1 and 2, respectively.
	bool **flag1 = InitializeFlag(xDim1, yDim1);
	bool **flag2 = InitializeFlag(xDim2, yDim2);
	//bool **unTouch = InitializeFlag(xDim1, yDim1);

	// Pre-computer the Confident Measure for two images.
	// Just for saving time
	bool** confidL = InitializeConfidL(im1, xDim1, yDim1, sThresh);
	bool** confidR = InitializeConfidL(im2, xDim2, yDim2, sThresh);

	// Initialize gradient confidence measures for two images.
	//double **grad1 = InitializeGrad(im1);
	//double **grad2 = InitializeGrad(im2);

	denseMatches.clear(); // Clear input dense matches.
	vector<Match> seed; // seed matches for propagation.
	seed.reserve(denseMatches.capacity());
	// list<Match> seed;
	vector<int> seedHeap;
	seedHeap.reserve(seed.capacity());

	// Initialize list of seed and dense matches by the input sparse matches.
	InitialHeap(seedHeap, matches.size());

	for (int i=0; i<matches.size(); ++i)
	{
		int x1 = (int)(matches[i].GetTarget().GetX() + 0.5);
		int y1 = (int)(matches[i].GetTarget().GetY() + 0.5);
		int x2 = (int)(matches[i].GetFound().GetX() + 0.5);
		int y2 = (int)(matches[i].GetFound().GetY() + 0.5);
		//double x1 = matches[i].GetTarget().GetX();
		//double y1 = matches[i].GetTarget().GetY();
		//double x2 = matches[i].GetFound().GetX();
		//double y2 = matches[i].GetFound().GetY();

		if (x1 < 0 || x1 > xDim1-1 || y1 < 0 || y1 > yDim1-1 || 
			x2 < 0 || x2 > xDim2-1 || y2 < 0 || y2 > yDim2-1)
		{
			continue;
		}

		// ZNCC correlation score.
		double score = matches[i].GetDist();
		if (score < -1.0 || score > 1.0)
		{
			continue;
		}

		Point p1(x1, y1);
		Point p2(x2, y2);

		Match m(p1, p2, score);
		m.SetIsKey(1); // Mark as an interest point.

		//seed.push_front(m);
		seed.push_back(m);
		denseMatches.push_back(m);
		
		// Update uniqueness flags.
		flag1[x1][y1] = false;
		flag2[x2][y2] = false;
	}

	// store candidate matching pair locally
	vector<Match> local;
	///////////////////////////////////////////////////////////////
	local.reserve(600);
	///////////////////////////////////////////////////////////////

	// Start propagation process.
	//while (seed.size() != 0)
	while(seedHeap.size() != 0)
	{
		int maxIdx = HeapMaxExtract(seed, seedHeap);
		Match s(seed[maxIdx]);
		local.clear();

		//vector<Match> local; // Store new candidate matches.

		int x1 = (int)(s.GetTarget().GetX());
		int y1 = (int)(s.GetTarget().GetY());
		int x2 = (int)(s.GetFound().GetX());
		int y2 = (int)(s.GetFound().GetY());

		// Propagate within the neighbor region of the seed.
		for (int i=-neighb; i<=neighb; ++i)
		{
			for (int j=-neighb; j<=neighb; ++j)
			{
				if (i == 0 && j == 0)
				{
					continue;
				}

				// each neighbor point in the first image.
				int u1 = x1 + i;
				int v1 = y1 + j;

				// center of the correspondence region in the second image.
				int uC = x2 + i;
				int vC = y2 + j;

				// exclude pixels out of boundary, 
				// 1 is for the gradient computation.
				if (u1 < 1 || u1 >= (xDim1-1) || 
					v1 < 1 || v1 >= (yDim1-1) || 
					uC < 1 || uC >= (xDim2-1) || 
					vC < 1 || vC >= (yDim2-1))
				{
					continue;
				}

				// Enforce the uniqueness constraint and 
				// gradient confidence measure of image 1.
				//if (!flag1[u1][v1] || 
				//	GradMeasure(u1, v1, im1) < sThresh)
				if (!flag1[u1][v1] || !confidL[u1][v1])// || !unTouch[u1][v1])
				{
					continue;
				}

				// corresponding point in the second image.
				int u2 = uC;
				int v2 = vC;
				double maxCorr = -1.0; // maximum ZNCC value.

				// Enforce the disparity gradient limit.
				for (int k=-dThresh; k<=dThresh; ++k)
				{
					for (int h=-dThresh; h<=dThresh; ++h)
					{
						int uN = uC + k;
						int vN = vC + h;

						// epipolar constraint
						if (!IsInLimit(u1, v1, uN, vN, fund, thresh))
						{
							continue;
						}

						//if (uN < Max(x2-neighb, 1) || 
						//	uN > Min(x2+neighb, xDim2-1) || 
						//	vN < Max(y2-neighb, 1) || 
						//	vN > Min(y2+neighb, yDim2-1))
						//{
						//	continue;
						//}
						if (uN < Maximun(x2-neighb, 1) || 
							uN > Minimun(x2+neighb, xDim2-1) || 
							vN < Maximun(y2-neighb, 1) || 
							vN > Minimun(y2+neighb, yDim2-1))
						{
							continue;
						}

						// modified by styx, Dec 17, 2007
						// for boundary problem--(here, radius = 2)
						//      so, x'+radius, x'-radius
						//          y'+radius, y'-radius 
						//      must be in the boundary of image2

						//if (uN < Max(x2-neighb, radius) || 
						//	uN > Min(x2+neighb, xDim2-1-radius) || 
						//	vN < Max(y2-neighb, radius) || 
						//	vN > Min(y2+neighb, yDim2-1-radius) ||
						//	!confidR[uN][vN] ||
						//	!flag2[uN][vN] )
						//{
						//	continue;
						//}


						// ZNCC value of (u1,v1) and (uN, vN).
						double corr = NormCrossCorr(u1, v1, 
							uN, vN, radius, im1, im2);
						
						if (corr > cThresh)
						{
							Point p1(u1, v1);
							Point p2(uN, vN);

							Match localTmp(p1, p2, corr);
							local.push_back(localTmp);
							// PushSeed2List(local, localTmp);
						}
					}
				} // end of for.
			}
		} // end of for.

		int localSize = local.size();
		
		// sort proagation by zncc score
		//sort(local.begin(), local.end(), less<Match>());

		// find the result
		sort(local.begin(), local.end(), less<Match>());
		for (int i = 0; i < localSize; i++)
		{
			//Match tmp(local.front());     // pull the ZNCC-best match(u, u') from(Local)
			//local.pop_front();
			Match tmp(local[i]);

			x1 = (int)tmp.GetTarget().GetX();
			y1 = (int)tmp.GetTarget().GetY();
			x2 = (int)tmp.GetFound().GetX();
			y2 = (int)tmp.GetFound().GetY();

			if (flag1[x1][y1] && flag2[x2][y2])
			{
				denseMatches.push_back(tmp);
					
				//seed.push_back(cand);
				//PushSeed2List(seed, tmp);
				seed.push_back(tmp);
				HeapInsert(seed, seedHeap, seed.size()-1);
				// Update uniqueness flags.
				flag1[x1][y1] = false;
				flag2[x2][y2] = false;
				//unTouch[x1][y1] = false;
			}
		}
	} // end of while.

	// Delete allocated uniqueness flags.
	if (flag1 != NULL)
	{
		for (int x=0; x<xDim1; x++)
		{
			delete [] flag1[x];
		}
		delete []flag1;
	}
	if (flag2 != NULL)
	{
		for (int x=0; x<xDim2; x++)
		{
			delete [] flag2[x];
		}
		delete []flag2;
	}

	//if (unTouch != NULL)
	//{
	//	for (int x=0; x<xDim1; x++)
	//	{
	//		delete [] unTouch[x];
	//	}
	//	delete []unTouch;
	//}
	
	// The confident measure matrix should be deleted
	if (confidL != NULL)
	{
		for (int x=0; x<xDim1; x++)
		{
			delete [] confidL[x];
		}
		delete []confidL;
	}
	if (confidR != NULL)
	{
		for (int x=0; x<xDim2; x++)
		{
			delete [] confidR[x];
		}
		delete []confidR;
	}

	/////
	//start = clock();
	// Sort final matches in a decreasing order of correlation coefficients.
	
	// sort(denseMatches.begin(), denseMatches.end(), less<Match>());
	
	/////
	//finish = clock();
	//time = (double)(finish - start) / CLOCKS_PER_SEC;
	//cout << "sorting uses " << time << "s" << endl;

	return denseMatches.size();
}



bool ** InitializeFlag(int xDim, int yDim)
{
	int x;
	// Allocate memory.
	bool **flag = new bool *[xDim];
	if (flag == NULL)
	{
		FatalError("InitializeFlag -- Allocating memory fails!");
	}
	for (x=0; x<xDim; ++x)
	{
		flag[x] = new bool[yDim];
		if (flag[x] == NULL)
		{
			FatalError("InitializeFlag -- Allocating memory fails!");
		}
	}

	// Initialize flags.
	for (x=0; x<xDim; ++x)
	{
		for (int y=0; y<yDim; ++y)
		{
			flag[x][y] = true;
		}
	}

	return flag;
}

bool** InitializeConfidL(Image& img, int xDim, int yDim, double sThresh)
{
	int x;
	// Allocate memory.
	bool **dgl = new bool* [xDim];
	if (dgl == 0)
	{
		FatalError("InitializeConfidL -- Allocating memory fails!");
	}
	for (x=0; x<xDim; ++x)
	{
		dgl[x] = new bool[yDim];
		if (dgl[x] == 0)
		{
			FatalError("InitializeConfidR -- Allocating memory fails!");
		}
	}

	// Initialize flags.
	for (x = 1; x < xDim-1; ++x)
	{
		for (int y = 1; y < yDim-1; ++y)
		{
			dgl[x][y] = (GradMeasure(x, y, img) > sThresh);
		}
	}

	return dgl;
}

/*
double ** InitializeGrad(Image &im)
{
	int xDim = im.GetXDim();
	int yDim = im.GetYDim();

	// Allocate memory.
	double **grad = new double *[xDim];
	if (grad == NULL)
	{
		FatalError("InitializeGrad -- Allocating memory fails!");
	}
	for (int x=0; x<xDim; ++x)
	{
		grad[x] = new double[yDim];
		if (grad[x] == NULL)
		{
			FatalError("InitializeGrad -- Allocating memory fails!");
		}
	}

	// Initialize gradient confidence measures.
	for (x=0; x<xDim; ++x)
	{
		for (int y=0; y<yDim; ++y)
		{
			if (x==0 || x==(xDim-1) || y==0 || y==(yDim-1))
			{
				grad[x][y] = -1.0;

				continue;
			}

			double pixel[4];
			pixel[0] = fabs(im(x, y-1) - im(x, y));
			pixel[1] = fabs(im(x, y+1) - im(x, y));
			pixel[2] = fabs(im(x-1, y) - im(x, y));
			pixel[3] = fabs(im(x+1, y) - im(x, y));

			double maxPixel = 0.0;
			for (int i=0; i<4; ++i)
			{
				if (pixel[i] > maxPixel)
				{
					maxPixel = pixel[i];
				}
			}

			grad[x][y] = maxPixel;
		}
	}

	return grad;
}
*/

/********************************************************************
   Compute confidence measure for each pixel x.
   s(x) = max{ |I(x+d)-I(x)|, d={(1,0),(-1,0),(0,1),(0,-1)} }.
********************************************************************/
double GradMeasure(int x, int y, Image &im)
{
	assert(x >= 1 && x < (im.GetXDim()-1) 
		&& y >= 1 && y < (im.GetYDim()-1));

	double pixel[4];
	pixel[0] = fabs(im(x, y-1) - im(x, y));
	pixel[1] = fabs(im(x, y+1) - im(x, y));
	pixel[2] = fabs(im(x-1, y) - im(x, y));
	pixel[3] = fabs(im(x+1, y) - im(x, y));

	double maxPixel = 0.0;
	for (int i=0; i<4; ++i)
	{
		if (pixel[i] > maxPixel)
		{
			maxPixel = pixel[i];
		}
	}

	return maxPixel;
}

/********************************************************************
   Resample input quasi-dense pixel correspondences to get 
   quasi-dense point correspondences with sub-pixel accuracy.
   Local geometric constraints are used, here we use affine transform.
   Resampled points include centers of patches and interest points.
********************************************************************/
int Resample(vector<Match> &sampledMatches, 
			 const vector<Match> &matches, 
			 Image &im1, Image &im2, int times,
			 vector<Match>* pMatch01, char* match12File)
{
	// sampledMatches means the common matches between image 0, 1 and 2
	// match12 store the matches between image 1 and 2 (but beside sampledMatches)
	vector<Match> match12;

	// store matches adjusted (matches between image0 and image1)
	vector<Match> matchNew01;

	int timetmp = (im1.GetXDim()/1200 + 1)*2;
	int square;
	if (timetmp == 2)
		square = timetmp * 4 - 2*times;
	else if (timetmp > 2)
		square = static_cast<int> (timetmp*4 / pow(2.0, times) + 0.1);

	// int square = (im1.GetXDim()/1200+1) * 8 - times*4; // divided square size.
	int reNum = 3;  // required number of matches for local fitting.
	//int radius = 2; // radius of ZNCC window.
	//double cThresh = 0.5; // ZNCC value threshold.

	int xDim1 = im1.GetXDim();	int yDim1 = im1.GetYDim();
	int xDim2 = im2.GetXDim();	int yDim2 = im2.GetYDim();

	// Divide image 1 into squares.
	int xSquare = xDim1 / square;       // square number of width.
	int ySquare = yDim1 / square;       // square number of height.
	int xMargin = (xDim1 % square) / 2; // margin of width after division.
	int yMargin = (yDim1 % square) / 2; // margin of height after division.

	int matchNum = matches.size();
	sampledMatches.clear();

	int squareNum = xSquare * ySquare;
	
	// local matches divided into squares.
	vector< vector<Match> > localMatches(squareNum);
	// local matches of interest points.
	vector< vector<Match> > localInterst(squareNum);
	vector< vector<Match> > match01Pat(squareNum);
	/////
	//int keyCount = 0;

	int i, j, n;
	// Assign all matches into local squares.
	AssignMatch2Local(localMatches, localInterst, matches, 
		xMargin, yMargin, xDim1, yDim1, xSquare, squareNum, square);

	// if pMatch01!=0, then assign matches between image0 and image1 into local patches
	if (pMatch01 != 0)
	{
		matchNew01.reserve(pMatch01->size());
		AssignMatch2Local(match01Pat, pMatch01, xMargin, yMargin, 
					   xDim1, yDim1, xSquare, squareNum, square);
	}

	//cout << keyCount << " keys." << endl;

	double xC1, yC1, xpre1, ypre1;
	int keyNum;
	bool flag;

	// Resample matches by local geometric constraints.
	for (j=0; j<ySquare; ++j)
	{
		for (i=0; i<xSquare; ++i)
		{
			int idx = i+j*xSquare;

			// local matches within current patch.
			vector<Match>& inMatches = localMatches[idx];
			// interest matches within current patch.
			vector<Match>& keyMatches = localInterst[idx];
			// pre-computed matches between image0 and image1
			vector<Match>& preMatches = match01Pat[idx];

			if (inMatches.size() <= reNum)
			{
				continue;
			}

			// Fit local affine transform in each patch.
			vector<int> inliers;
			////////////////////////////////////////////////////
			inliers.reserve(inMatches.size());
			///////////////////////////////////////////////////

			NyxMat<double> affine;
			double thresh = 0.001;

			int flag = RansacFitAffine(affine, inliers, inMatches, thresh, 30, 30);
			if (flag == 0)
			{
				//cout << "RANSAC failed." << endl;
				continue;
			}

			// store the center of the patch
			xC1 = xMargin + i*square + square/2.0 - 0.5;
			yC1 = yMargin + j*square + square/2.0 - 0.5;

			if (pMatch01 == 0)	// use local center of current patch in image 1.
			{
				AddPnt2Vec(xC1, yC1, im1, im2, false, affine, sampledMatches);
			}
			else				//-- find common matches between 3 images.		  
			{
				keyNum = preMatches.size();
				if (keyNum == 0) // no pre-computed matches, then use the center of the patch
				{
					AddPnt2Vec(xC1, yC1, im1, im2, false, affine, match12);
				}
				else
				{
					for (n=0; n<keyNum; ++n)
					{
						// interest point within current patch in image 1.
						xpre1 = preMatches[n].GetFound().GetX();
						ypre1 = preMatches[n].GetFound().GetY();

						// Compute corresponding interest point in image 2.
						flag = AddPnt2Vec(xpre1, ypre1, im1, im2, false, affine, sampledMatches);
						if (flag)
						{
							matchNew01.push_back(preMatches[n]);
						}

					} // end of interest matches.
				}
			}

			// Add all corresponding points of interest within the patch.
			keyNum = keyMatches.size();
			for (n=0; n<keyNum; ++n)
			{
				// interest point within current patch in image 1.
				xpre1 = keyMatches[n].GetTarget().GetX();
				ypre1 = keyMatches[n].GetTarget().GetY();

				// Compute corresponding interest point in image 2.
				if (pMatch01 == 0)
				{
					flag = IsSamePoint(xpre1, ypre1, xC1, yC1);
					if (flag == false)   // key point is not same as center point
						// AddPnt2Vec(xpre1, ypre1, im1, im2, true, affine, sampledMatches);
						sampledMatches.push_back(keyMatches[n]); // insterest keys are not mapped
				}
				else
				{
					if (preMatches.size() == 0)
					{
						flag = IsSamePoint(xpre1, ypre1, xC1, yC1);
						if (flag == false)
							// AddPnt2Vec(xpre1, ypre1, im1, im2, true, affine, match12);
							sampledMatches.push_back(keyMatches[n]);
					}
					else
					{
						// judge whether the key interest is between matches of image0 and image1
						flag = IsInMatchVec(xpre1, ypre1, preMatches);
						if (flag == false)
							// AddPnt2Vec(xpre1, ypre1, im1, im2, true, affine, match12);
							sampledMatches.push_back(keyMatches[n]); 
					}
				}
			} // end of interest matches.
		}
	} // end of patches.

	// write out the matches between image 1 and image 2 and at same time not the matches
	// between image 0 and image 1
	if (pMatch01 != 0)
	{
		(*pMatch01) = matchNew01;   // copy the adjusted matches (sorted and partially deleted)
			
		assert(match12File);
		WriteMatchFile(match12File, match12);  		 //corruption?
	}
	
//	printf("Corruption???????\n");

	return sampledMatches.size();
}



// find corresponding points using affine matrix, and then store that match into a vector
bool AddPnt2Vec(double x1, double y1, Image &im1, Image &im2,
				bool isKey, NyxMat<double> &affine, vector<Match>& matches)
{
	int xDim2 = im2.GetXDim();
	int yDim2 = im2.GetYDim();

	// compute corresponding point
	double x2 = 
		affine(0,0)*x1 + affine(0,1)*y1 + affine(0,2);
	double y2 = 
		affine(1,0)*x1 + affine(1,1)*y1 + affine(1,2);

	if (x2 >= 0 && x2 <  double(xDim2)-0.5 && y2 >= 0 && y2 < double(yDim2)-0.5)
	{
		Point p1(x1, y1);
		Point p2(x2, y2);

		// ZNCC value of (x1,y1) and (x2, y2).
		double corr = NormCrossCorr(
			(int)(x1+0.5), (int)(y1+0.5), 
			(int)(x2+0.5), (int)(y2+0.5), 
			PROPA_RADIUS, im1, im2);

		if (corr > C_THRESH)
		{
			Match m(p1, p2, corr);
			m.SetIsKey(isKey); // Mark as an interest point.

			matches.push_back(m);
			return true;
		}
	}
	return false;
}


// judge whether two points are same
bool IsSamePoint(Point& pnt1, Point& pnt2)
{
	if (fabs(pnt1.GetX()-pnt2.GetX()) < DIST_SAME_THRESH && 
		fabs(pnt1.GetY()-pnt2.GetY()) < DIST_SAME_THRESH)
	{
		return true;
	}

	return false;
}

// overload 
bool IsSamePoint(double x1, double y1, double x2, double y2)
{
	if (fabs(x1-x2) < DIST_SAME_THRESH && 
		fabs(y1-y2) < DIST_SAME_THRESH)
	{
		return true;
	}

	return false;
}


// judge whether pnt1 is included in vector "matches"
bool IsInMatchVec(Point& pnt1, vector<Match>& matches)
{
	int featureNum = matches.size();
	Point secp;

	for (int i = 0; i < featureNum; ++i)
	{
	  secp = matches[i].GetFound();
		if (IsSamePoint(pnt1, secp))
		{
			return true;
		}
	}

	return false;
}

// overload 
bool IsInMatchVec(int x1, int x2, vector<Match>& matches)
{
	int featureNum = matches.size();
	
	for (int i = 0; i < featureNum; ++i)
	{
		if (IsSamePoint(x1, x2, matches[i].GetFound().GetX(), matches[i].GetFound().GetY()))
		{
			return true;
		}
	}

	return false;
}


// assign matches into match queues of local patches
void AssignMatch2Local(vector<vector<Match> >& localMatches,
					   vector<vector<Match> >& localInterst,
					   const vector<Match>& matches, int xMargin, int yMargin, 
					   int xDim1, int yDim1, int xSquare, int squareNum, int square)
{
	int n;
	int matchNum = matches.size();
	int x1, y1;

	// Assign all matches into local squares.
	for (n=0; n<matchNum; ++n)
	{
		x1 = (int)matches[n].GetTarget().GetX();
		y1 = (int)matches[n].GetTarget().GetY();

		// Get rid of points in the margins.
		if (x1 < xMargin || x1 >= (xDim1-xMargin) || 
			y1 < yMargin || y1 >= (yDim1-yMargin))
		{
			continue;
		}

		// local square index.
		int idx = (int)((x1 - xMargin) / square) + 
			xSquare * (int)((y1 - yMargin) / square);

		if (idx >= 0 && idx < squareNum)
		{
			localMatches[idx].push_back(matches[n]);

			// Store interest matches.
			if (matches[n].GetIsKey() == 1)
			{
				localInterst[idx].push_back(matches[n]);

				/////
				//++keyCount;
			}
		}
	}
}


// overload
void AssignMatch2Local(vector<vector<Match> >& match01Pat,
					   vector<Match>* pMatch01, int xMargin, int yMargin, 
					   int xDim1, int yDim1, int xSquare, int squareNum, int square)
{
	int n;
	int matchNum = pMatch01->size();
	int x1, y1;

	// match01Pat.reserve(squareNum);
	for (n=0; n<matchNum; ++n)
	{
		x1 = (int)(*pMatch01)[n].GetFound().GetX();
		y1 = (int)(*pMatch01)[n].GetFound().GetY();

		// Get rid of points in the margins.
		if (x1 < xMargin || x1 >= (xDim1-xMargin) || 
			y1 < yMargin || y1 >= (yDim1-yMargin))
		{
			continue;
		}

		// local square index.
		int idx = (int)((x1 - xMargin) / square) + 
			xSquare * (int)((y1 - yMargin) / square);

		if (idx >= 0 && idx < squareNum)
		{
			match01Pat[idx].push_back((*pMatch01)[n]);
		}
	}
}
