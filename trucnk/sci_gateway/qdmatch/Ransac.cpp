//////////////////////////////////////////////////////////////////////////////////////////////////
// Author:		Styx(Zhenhui Xu)
// Version:		0.1
// Date:		May 20, 2009 
// Description: This file is based on the Matlab routines by Peter Kovesi. 
// Others:      available at http://www.csse.uwa.edu.au/~pk/Research/MatlabFns/.
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


#include <ctime>
#include "Ransac.h"
#include "NyxSvd.h"


/********************************************************************
   'm' -- the model having the greatest number of inliers.
   'inliers' -- an array of indices of the elements in the data that  
                are inliers for the best model.
   'points1' -- input matched points from the first image.
   'points2' -- input matched points from the second image.
   'n' -- the minimum number of samples from the data required to 
          fit a model.
   't' -- the distance threshold between a data point and the model 
          used to decide whether the point is an inlier or not.
   'fit' -- the fitting model, can be 'eFundmatrix' or 'eHomography'.
   'maxTrials' -- maximum number of trials before we give up.
   'maxDataTrials' -- max number of attempts to select a non-degenerate data set.

   This function robustly fits a model to data with the RANSAC algorithm.
********************************************************************/
int Ransac(NyxMat<double> &m, vector<int> &inliers, 
		   const NyxMat<double> &points1, const NyxMat<double> &points2, 
		   int n, double t, FitModel fit, int maxTrials, int maxDataTrials)
{
	int dim = points1.GetRow();
	int matchNum = points1.GetCol();
	
	double p = 0.95;         // desired probability of choosing at least 
	                         // one sample free from outliers.

	NyxMat<double> bestM;    // best fitted model matrix.
	vector<int> bestInliers; // best inliers.
	int bestScore = 0;       // best inliers number.
	int trialCount = 0;      // count for trails.
	int k = 1;               // initial number of trials.

	int i, j;
	int count;
	bool degenerate;

	// reset the seed of random generator
	srand((unsigned)time(0)); 
	while (k > trialCount)
	{
		count = 1;
		degenerate = true;
		
		while (degenerate)
		{
			// Select random n points to form a trial model.
			NyxMat<int> ind(1, n, 0);

			// create the index --> [0 : matchNum-1]
			NyxMat<int> index(1, matchNum, GRIDY);
			
			// reset the seed of random generator
			// srand((unsigned)time(0)); 
		    for (i = 0; i < n; ++i)
			{
				// Generate a random integer belong to [0, matchNum-i)
				int aRand = rand() % (matchNum - i);

				ind(0, i) = index(0, aRand);
				index(0, aRand) = index(0, matchNum - 1 - i);
			}

			NyxMat<double> drawnPoints1(dim, n, 0);
			NyxMat<double> drawnPoints2(dim, n, 0);

			for (j = 0; j < n; ++j)
			{
				for (i = 0; i < dim; ++i)
				{
					int jj = ind(0, j);
					drawnPoints1(i, j) = points1(i, jj);
					drawnPoints2(i, j) = points2(i, jj);
				}
			}
			
			// Check whether these points are degenerate.
			degenerate = DegenFun(drawnPoints1, drawnPoints2, fit);

			// Fit model to this random selection of data points.
			if (!degenerate)
			{
				m = FittingFun(drawnPoints1, drawnPoints2, fit);
				
				if (m.GetMatrix() == 0)
				{
					degenerate = true;
				}
			}

			++count;

			// Too many trials for selecting random points.
			if (count > maxDataTrials)
			{
				// cerr << "unable to select a nondegenerate data set." << endl;
				return 0;
			}
		}

		// Evaluate distances between points and model.
		//inliers.clear();
		DistFun(inliers, m, points1, points2, t, fit);

		// Record the largest set of inliers so far.
		int inlierNum = inliers.size();
		if (inlierNum > bestScore)
		{
			bestScore = inlierNum;
			bestInliers = inliers;
			bestM = m;

			// Update estimation of the needed number of trials.
			double inlierFrac = (double)inlierNum / (double)matchNum;
			double pNoOutliers = 1 - pow(inlierFrac, n);
			//pNoOutliers = Max(EPS, pNoOutliers);
			//pNoOutliers = Min(1-EPS, pNoOutliers);
			pNoOutliers = Maximun(EPS, pNoOutliers);
			pNoOutliers = Minimun(1-EPS, pNoOutliers);
			k = (int)(log(1-p) / log(pNoOutliers)) + 1;
		}
		
		++trialCount;

		// Too many trials.
		if (trialCount > maxTrials)
		{
			//cout << "ransac reached the maximum number of " 
			//	<< maxTrials << " trials." << endl;
			break;
		}
	}

	// Return the best model and inliers.
	if (bestM.GetMatrix() != 0)
	{
		m = bestM;
		inliers = bestInliers;
		return 1;
	}
	else
	{
		// cout << "ransac was unable to find a useful solution." << endl;
		return 0;
	}
}


/********************************************************************
   This function evaluates the distances from the model 'm' to data 
   'points1' and 'points2'.
   Return the indices of elements in the data that are inliers, 
   that is, the points that are within distance 't' of the model.
********************************************************************/
void DistFun(vector<int> &inliers, const NyxMat<double> &m, 
			 const NyxMat<double> &points1, const NyxMat<double> &points2, 
			 double t, FitModel fit)
{
	assert(fit == eFundmatrix || fit == eHomography || fit == eAffine);

	if (fit == eFundmatrix)
	{
		FundDist(inliers, m, points1, points2, t);
	} 
	else if (fit == eHomography)
	{
		HomogDist(inliers, m, points1, points2, t);
	}
	else
	{
		AffineDist(inliers, m, points1, points2, t);
	}
}

/********************************************************************
   This function fits a model to data 'points1' and 'points2'.
********************************************************************/
NyxMat<double> FittingFun(const NyxMat<double> &points1, 
						  const NyxMat<double> &points2, FitModel fit)
{
	assert(fit == eFundmatrix || fit == eHomography || fit == eAffine);

	if (fit == eFundmatrix)
	{
		return Fundmatrix(points1, points2);
	}
	else if (fit == eHomography)
	{
		return Homography(points1, points2);
	}
	else
	{
		return Affine(points1, points2);
	}
}

/********************************************************************
   This function determines whether a set of data points will 
   produce a degenerate model.
   This is used to discard random samples that do not result in 
   useful models.
   If we cannot determine the degeneracy explicitly, we can just 
   return 'false'.
********************************************************************/
bool DegenFun(const NyxMat<double> &points1, const NyxMat<double> &points2, FitModel fit)
{
	assert(fit == eFundmatrix || fit == eHomography || fit == eAffine);

	if (fit == eFundmatrix)
	{
		return false;
	}
	else if (fit == eHomography)
	{
		return IsDegenHomog(points1, points2);
	}
	else
	{
		return IsDegenAffine(points1, points2);
	}
}



/********************************************************************
   'model' -- the 3x3 fundamental matrix such that P2'*F*P1 = 0.
   'inliers' -- an array of indices of the elements in the data that  
                are inliers for the best model.
   'matches' -- input matched points from two images.
   't' -- The distance threshold between a data point and the model 
          used to decide whether the point is an inlier or not.
   This function fits fundamental matrix using RANSAC.
********************************************************************/
int RansacFitFund(NyxMat<double> &model, vector<int> &inliers, 
				  const vector<Match> &matches, double t)
{
	int i;

	// clear model
	model = NyxMat<double>(3, 3, 0);
	//for (int i = 0; i < model.GetRow(); ++i)
	//	for (int j = 0; j < model.GetCol(); ++j)
	//		model(i, j) = 0;

	inliers.clear();

	int matchNum = matches.size();

	if (matchNum < 8) // can't find an fundamental matrix, return error
		return 0;
	
	NyxMat<double> points1(3, matchNum, 0);
	NyxMat<double> points2(3, matchNum, 0);
	
	for (i = 0; i < matchNum; ++i)
	{
		points1(0, i) = matches[i].GetTarget().GetX();
		points1(1, i) = matches[i].GetTarget().GetY();
		points1(2, i) = 1;
		
		points2(0, i) = matches[i].GetFound().GetX();
		points2(1, i) = matches[i].GetFound().GetY();
		points2(2, i) = 1;
	}

	// Normalize each set of points so that the origin is at centroid 
	// and mean distance from origin is sqrt(2). 
	NyxMat<double> trans1, trans2;
	points1 = NormalizePts(trans1, points1);
	points2 = NormalizePts(trans2, points2);

	// The minimum number of samples required to fit a fundamental 
	// matrix model is 7, here we use 8-points solution.
	int n = 8;
	NyxMat<double> m(3, 3, 0);

	int flag = Ransac(m, inliers, points1, points2, n, t, eFundmatrix);

	if (flag == 0)
	{
		//cout << "ransac was unable to find a useful solution." << endl;
		return flag;
	}

	int inlierNum = inliers.size();
	NyxMat<double> inlierPts1(3, inlierNum, 0);
	NyxMat<double> inlierPts2(3, inlierNum, 0);

	for (i=0; i<inlierNum; ++i)
	{
		inlierPts1(0, i) = points1(0, inliers[i]);
		inlierPts1(1, i) = points1(1, inliers[i]);
		inlierPts1(2, i) = points1(2, inliers[i]);

		inlierPts2(0, i) = points2(0, inliers[i]);
		inlierPts2(1, i) = points2(1, inliers[i]);
		inlierPts2(2, i) = points2(2, inliers[i]);
	}

	m = Fundmatrix(inlierPts1, inlierPts2);

	// Denormalize : F = T2'*F*T1
	model = (!trans2)*m*trans1;

	return flag;
}

/********************************************************************
   Compute fundamental matrix from 8 or more points.
********************************************************************/
NyxMat<double> Fundmatrix(const NyxMat<double> &points1, const NyxMat<double> &points2)
{
	// Normalize each set of points.
	NyxMat<double> trans1, trans2;
	NyxMat<double> newPts1 = NormalizePts(trans1, points1);
	NyxMat<double> newPts2 = NormalizePts(trans2, points2);

	int i;
	// Build the constrain matrix.
	int ptsNum = points1.GetCol();

	NyxMat<double> a(ptsNum, 9, 0);

	for (i = 0; i < ptsNum; ++i)
	{
		a(i, 0) = newPts2(0, i) * newPts1(0, i);
		a(i, 1) = newPts2(0, i) * newPts1(1, i);
		a(i, 2) = newPts2(0, i);
		a(i, 3) = newPts2(1, i) * newPts1(0, i);
		a(i, 4) = newPts2(1, i) * newPts1(1, i);
		a(i, 5) = newPts2(1, i);
		a(i, 6) = newPts1(0, i);
		a(i, 7) = newPts1(1, i);
		a(i, 8) = 1;
	}

	//Mm transAA = mtimes(ctranspose(a), a);
	//Mm u, d, v;
	//eig(transAA, i_o, v, d);

	// SVD decomposing
	NyxMat<double> u, d, v;
	nyx_svd(a, 0, 1, u, d, v);
	
	// Extract fundamental matrix from the column of V 
	// corresponding to smallest singular value.
	NyxMat<double> fundM(3, 3, 0);
	for (i=0; i<3; ++i)
	{
		for (int j=0; j<3; ++j)
		{
			fundM(i, j) = v(i*3+j, 8);
		}
	}

	// Enforce constraint that fundamental matrix has rank 2 by performing 
	// a SVD and then reconstructing with the two largest singular values.
	nyx_svd(fundM, 1, 1, u, d, v);
	NyxMat<double> diagD(3, 3, 0);

	// cout << d;

	diagD(0, 0) = d(0, 0);
	diagD(1, 1) = d(1, 0);
	fundM = u*diagD*(!v);

	// Denormalize : F = T2'*F*T1.
	fundM = (!trans2)*fundM*trans1;

	return fundM;
}

/********************************************************************
   Evaluate the first order approximation of the geometric error 
   (Sampson distance) of the fit of a fundamental matrix with 
   respect to a set of matched points as needed by RANSAC.
********************************************************************/
void FundDist(vector<int> &inliers, const NyxMat<double> &fundM, 
			  const NyxMat<double> &points1, const NyxMat<double> &points2, double t)
{
	inliers.clear();

	int ptsNum = points1.GetCol();

	NyxMat<double> p2tFp1(1, ptsNum, 0);
	NyxMat<double> tmp;

	for (int i = 0; i < ptsNum; ++i)
	{
		NyxMat<double> pli(3, 1, 0);
		pli(0, 0) = points1(0, i);
		pli(1, 0) = points1(1, i);
		pli(2, 0) = points1(2, i);

		NyxMat<double> p2i(1, 3, 0);
		p2i(0, 0) = points2(0, i);
		p2i(0, 1) = points2(1, i);
		p2i(0, 2) = points2(2, i);

		tmp = p2i*fundM*pli;
		p2tFp1(0, i) = tmp(0,0);
	}

	NyxMat<double> fP1(fundM * points1);	   // F * P1
	NyxMat<double> ftP2((!fundM) * points2);   // F' * P2

	// Evaluate distances.
	for (int i = 0; i < ptsNum; ++i)
	{
		double dist = p2tFp1(0, i) * p2tFp1(0, i);
		dist /= fP1(0, i) * fP1(0, i) + fP1(1, i) * fP1(1, i) + 
			ftP2(0, i) * ftP2(0, i) + ftP2(1, i) * ftP2(1, i);

		if (fabs(dist) < t)
		{
			inliers.push_back(i);
		}
	}
}


/*--------------------------------Normalization------------------------------*/

/********************************************************************
   This function translates and normalizes a set of 2D homogeneous 
   points so that their centroid is at the origin and their 
   mean distance from the origin is sqrt(2).
   This process typically improves the conditioning of any equations 
   used to solve homographies, fundamental matrices etc.
********************************************************************/
NyxMat<double> NormalizePts(NyxMat<double> &trans, const NyxMat<double> &pts)
{
	int dim = pts.GetRow();
	int num = pts.GetCol();

	if (dim != 3)
		FatalError("pts must be 3xn.");

	NyxMat<double> newPts(pts);
	int i;
	for (i = 0; i < num; ++i)
	{
		assert(pts(2,i) != 0); // infinite point.

		newPts(0, i) = newPts(0, i) / newPts(2, i);
		newPts(1, i) = newPts(1, i) / newPts(2, i);
		newPts(2, i) = 1;
	}

	double center[2]={0, 0}; // centroid of finite points.

	for (i = 0; i < num; ++i)
	{
		center[0] += newPts(0, i);
		center[1] += newPts(1, i);
	}
	center[0] /= num;
	center[1] /= num;

	double meanDist = 0.0;
	
	for (i = 0; i < num; ++i)
	{
		meanDist += sqrt((newPts(0, i)-center[0])*(newPts(0, i)-center[0]) + 
			(newPts(1, i)-center[1])*(newPts(1, i)-center[1]));
	}
	meanDist /= num;

	double scale = sqrt(static_cast<double>(2)) / meanDist;

	trans = NyxMat<double>(3, 3, 0);
	trans(0, 0) = scale;
	trans(0, 2) = - scale * center[0];
	trans(1, 1) = scale;
	trans(1, 2) = - scale * center[1];
	trans(2, 2) = 1;

	return trans * newPts;
}


/********************************************************************
   Normalize array of homogeneous coordinates to a scale of 1.
********************************************************************/
NyxMat<double> NormalizeCrd(const NyxMat<double> &pts)
{
	int dim = pts.GetRow();
	int num = pts.GetCol();
	
	NyxMat<double> newPts(pts);

	for (int i = 0; i < num; ++i)
	{
		if (fabs(pts(dim-1, i)) < EPS)
		{
			cout << "some point is at infinity." << endl;
			continue;
		}

		for (int j = 0; j < dim-1; ++j)
		{
			newPts(j, i) = pts(j, i) / pts(dim-1, i);
		}
		newPts(dim-1, i) = 1;
	}

	return newPts;
}


/*-----------------------RANSAC fitting homography model---------------------*/

/********************************************************************
   'model' -- the 3x3 homography such that P2 = H*P1.
   'inliers' -- an array of indices of the elements in the data that  
                are inliers for the best model.
   'matches' -- input matched points from two images.
   't' -- The distance threshold between a data point and the model 
          used to decide whether the point is an inlier or not.
   This function fits homography using RANSAC.
********************************************************************/
int RansacFitHomog(NyxMat<double> &model, vector<int> &inliers, 
				   const vector<Match> &matches, double t)
{
	int i;
	model = NyxMat<double>(3, 3, 0);

	inliers.clear();

	int matchNum = matches.size();
	NyxMat<double> points1(3, matchNum, 0);
	NyxMat<double> points2(3, matchNum, 0);
	
	for (i = 0; i < matchNum; ++i)
	{
		points1(0, i) = matches[i].GetTarget().GetX();
		points1(1, i) = matches[i].GetTarget().GetY();
		points1(2, i) = 1;
		
		points2(0, i) = matches[i].GetFound().GetX();
		points2(1, i) = matches[i].GetFound().GetY();
		points2(2, i) = 1;
	}

	// Normalize each set of points.
	NyxMat<double> trans1, trans2;
	points1 = NormalizePts(trans1, points1);
	points2 = NormalizePts(trans2, points2);

	// Minimum number of points needed to fit a homography is 4.
	int n = 4;
	NyxMat<double> m(3, 3, 0);

	int flag = Ransac(m, inliers, points1, points2, n, t, eHomography);

	if (flag == 0)
	{
		// std::cout << "ransac was unable to find a useful solution." << endl;
		return flag;
	}

	int inlierNum = inliers.size();
	NyxMat<double> inlierPts1(3, inlierNum, 0);
	NyxMat<double> inlierPts2(3, inlierNum, 0);

	for (i = 0; i < inlierNum; ++i)
	{
		inlierPts1(0, i) = points1(0, inliers[i]);
		inlierPts1(1, i) = points1(1, inliers[i]);
		inlierPts1(2, i) = points1(2, inliers[i]);

		inlierPts2(0, i) = points2(0, inliers[i]);
		inlierPts2(1, i) = points2(1, inliers[i]);
		inlierPts2(2, i) = points2(2, inliers[i]);
	}

	// Now do a final least squares fit on the 
	// data points considered to be inliers.
	m = Homography(inlierPts1, inlierPts2);

	// Denormalize : H = T2\H*T1
	model = MatInv((!trans2)*trans2)*(!trans2)*m*trans1;

	return flag;
}

/********************************************************************
   Evaluate the symmetric transfer error of a homography with 
   respect to a set of matched points as needed by RANSAC.
********************************************************************/
void HomogDist(vector<int> &inliers, const NyxMat<double> &homog, 
			   const NyxMat<double> &points1, const NyxMat<double> &points2, double t)
{
	inliers.clear();

	// Calculate, in both directions, the transfered points.
	NyxMat<double> hP1(homog * points1);			// H*P1
	NyxMat<double> invhP2(MatInv((!homog)*homog)*(!homog)*points2); // H^(-1)*P2 or H\P2

	// Normalize: make the homogeneous scales for all points are 1.
	NyxMat<double> newPts1 = NormalizeCrd(points1);
	NyxMat<double> newPts2 = NormalizeCrd(points2);
	
	hP1 = NormalizeCrd(hP1);
	invhP2 = NormalizeCrd(invhP2);

	int ptsNum = points1.GetCol();

	for (int i = 0; i < ptsNum; ++i)
	{
		double dist = 
			pow((newPts1(0, i) - invhP2(0, i)), 2) + 
			pow((newPts1(1, i) - invhP2(1, i)), 2) + 
			pow((newPts2(0, i) -    hP1(0, i)), 2) + 
			pow((newPts2(1, i) -    hP1(1, i)), 2);

		if (fabs(dist) < t)
		{
			inliers.push_back(i);
		}
	}
}


/********************************************************************
   Compute homography from 4 or more points.
   points1, points2: 3xN
********************************************************************/
NyxMat<double> Homography(const NyxMat<double>& points1, const NyxMat<double>& points2)
{
	// Normalize each set of points.
	int dim = points1.GetRow();
	int ptsNum = points1.GetCol();
	NyxMat<double> trans1, trans2;
	NyxMat<double> newPts1 = NormalizePts(trans1, points1);
	NyxMat<double> newPts2 = NormalizePts(trans2, points2);

	NyxMat<double> a(3*ptsNum, 9, 0);

	for (int i=0; i<ptsNum; ++i)
	{
		a(i*3, 3) = -newPts1(0, i) * newPts2(2, i);
		a(i*3, 4) = -newPts1(1, i) * newPts2(2, i);
		a(i*3, 5) = -newPts1(2, i) * newPts2(2, i);

		a(i*3, 6) = newPts1(0, i) * newPts2(1, i);
		a(i*3, 7) = newPts1(1, i) * newPts2(1, i);
		a(i*3, 8) = newPts1(2, i) * newPts2(1, i);

		a(i*3+1, 0) = newPts1(0, i) * newPts2(2, i);
		a(i*3+1, 1) = newPts1(1, i) * newPts2(2, i);
		a(i*3+1, 2) = newPts1(2, i) * newPts2(2, i);

		a(i*3+1, 6) = -newPts1(0, i) * newPts2(0, i);
		a(i*3+1, 7) = -newPts1(1, i) * newPts2(0, i);
		a(i*3+1, 8) = -newPts1(2, i) * newPts2(0, i);

		a(i*3+2, 0) = -newPts1(0, i) * newPts2(1, i);
		a(i*3+2, 1) = -newPts1(1, i) * newPts2(1, i);
		a(i*3+2, 2) = -newPts1(2, i) * newPts2(1, i);

		a(i*3+2, 3) = newPts1(0, i) * newPts2(0, i);
		a(i*3+2, 4) = newPts1(1, i) * newPts2(0, i);
		a(i*3+2, 5) = newPts1(2, i) * newPts2(0, i);
	}


	//Mm transAA = mtimes(ctranspose(a), a);
	//Mm d, v;
	//eig(transAA, i_o, v, d);

	NyxMat<double> u, d, v;
	nyx_svd(a, 0, 1, u, d, v);

	// Extract homography
	NyxMat<double> homog(3, 3, 0);
	for (int i=0; i<3; ++i)
	{
		for (int j=0; j<3; ++j)
		{
			homog(i, j) = v(i*3+j, 8);
		}
	}

	// Denormalize : H = T2\H*T1
	homog = MatInv((!trans2)*trans2)*(!trans2)*homog*trans1;

	return homog;
}

/********************************************************************
   Determine if a set of 4 pairs of matched points give rise to a 
   degeneracy in the calculation of a homography as needed by RANSAC.
   This involves testing whether any 3 of the 4 points in each set 
   is collinear. 
********************************************************************/
bool IsDegenHomog(const NyxMat<double>& points1, const NyxMat<double> &points2)
{
	assert(points1.GetRow() == 3 && points2.GetRow() == 3);
	assert(points1.GetCol() == 4 && points2.GetCol() == 4);

	NyxMat<double> p11(1, 3, 0);
	p11(0, 0) = points1(0, 0);
	p11(0, 1) = points1(1, 0);
	p11(0, 2) = points1(2, 0);

	NyxMat<double> p12(1, 3, 0);
	p12(0, 0) = points1(0, 1);
	p12(0, 1) = points1(1, 1);
	p12(0, 2) = points1(2, 1);

	NyxMat<double> p13(1, 3, 0);
	p13(0, 0) = points1(0, 2);
	p13(0, 1) = points1(1, 2);
	p13(0, 2) = points1(2, 2);

	NyxMat<double> p14(1, 3, 0);
	p14(0, 0) = points1(0, 3);
	p14(0, 1) = points1(1, 3);
	p14(0, 2) = points1(2, 3);

	NyxMat<double> p21(1, 3, 0);
	p21(0, 0) = points2(0, 0);
	p21(0, 1) = points2(1, 0);
	p21(0, 2) = points2(2, 0);

	NyxMat<double> p22(1, 3, 0);
	p22(0, 0) = points2(0, 1);
	p22(0, 1) = points2(1, 1);
	p22(0, 2) = points2(2, 1);

	NyxMat<double> p23(1, 3, 0);
	p23(0, 0) = points2(0, 2);
	p23(0, 1) = points2(1, 2);
	p23(0, 2) = points2(2, 2);

	NyxMat<double> p24(1, 3, 0);
	p24(0, 0) = points2(0, 3);
	p24(0, 1) = points2(1, 3);
	p24(0, 2) = points2(2, 3);
	
	bool isDeg = 
		IsCollinear(p11, p12, p13) || 
		IsCollinear(p11, p12, p14) || 
		IsCollinear(p11, p13, p14) || 
		IsCollinear(p12, p13, p14) || 
		IsCollinear(p21, p22, p23) || 
		IsCollinear(p21, p22, p24) || 
		IsCollinear(p21, p23, p24) || 
		IsCollinear(p22, p23, p24);

	return isDeg;
}

/********************************************************************
   Check if the input three points are collinear.
********************************************************************/
bool IsCollinear(const NyxMat<double> &p1, 
				 const NyxMat<double> &p2, 
				 const NyxMat<double> &p3)
{
	if (Mat2Norm(CrossMulti(p2-p1, p3-p1)) < EPS)
	{
		return true;
	}
	else
	{
		return false;
	}
}

/*-------------------------RANSAC fitting affine model-----------------------*/

/********************************************************************
   'model' -- the 2x3 affine transform such that P2 = A*P1.
   'inliers' -- an array of indices of the elements in the data that  
                are inliers for the best model.
   'matches' -- input matched points from two images.
   't' -- The distance threshold between a data point and the model 
          used to decide whether the point is an inlier or not.
   This function fits affine transform using RANSAC.
********************************************************************/
int RansacFitAffine(NyxMat<double> &model, vector<int> &inliers, 
					const vector<Match> &matches, double t, 
					int maxTrials, int maxDataTrials)
{
	int i;
	model = NyxMat<double>(2, 3, 0);

	inliers.clear();

	int matchNum = matches.size();
	NyxMat<double> points1(3, matchNum, 0);
	NyxMat<double> points2(3, matchNum, 0);
	
	for (i=0; i<matchNum; ++i)
	{
		points1(0, i) = matches[i].GetTarget().GetX();
		points1(1, i) = matches[i].GetTarget().GetY();
		points1(2, i) = 1;
		
		points2(0, i) = matches[i].GetFound().GetX();
		points2(1, i) = matches[i].GetFound().GetY();
		points2(2, i) = 1;
	}

	// Minimum number of points needed to fit an affine transform is 3.
	int n = 3;
	NyxMat<double> m(2, 3, 0);

	int flag = Ransac(m, inliers, points1, points2, n, t, eAffine, maxTrials, maxDataTrials);

	if (flag == 0)
	{
		// std::cout << "ransac was unable to find a useful solution." << endl;
		return flag;
	}

	int inlierNum = inliers.size();
	NyxMat<double> inlierPts1(3, inlierNum, 0);
	NyxMat<double> inlierPts2(3, inlierNum, 0);

	for (i=0; i<inlierNum; ++i)
	{
		inlierPts1(0, i) = points1(0, inliers[i]);
		inlierPts1(1, i) = points1(1, inliers[i]);
		inlierPts1(2, i) = points1(2, inliers[i]);

		inlierPts2(0, i) = points2(0, inliers[i]);
		inlierPts2(1, i) = points2(1, inliers[i]);
		inlierPts2(2, i) = points2(2, inliers[i]);
	}

	// Now do a final least squares fit on the 
	// data points considered to be inliers.
	model = Affine(inlierPts1, inlierPts2);

	return flag;
}



/********************************************************************
   Evaluate the symmetric transfer error of an affine transform with 
   respect to a set of matched points as needed by RANSAC.
   Here, we use the same as homography, whose value is set to :
     a11 a12 a13
     a21 a22 a23
     0   0   1
********************************************************************/
void AffineDist(vector<int> &inliers, const NyxMat<double> &affine, 
				const NyxMat<double> &points1, const NyxMat<double> &points2, double t)
{
	inliers.clear();

	NyxMat<double> homog(3, 3, 0);
	for (int i = 0; i < 2; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			homog(i, j) = affine(i,j);
		}
	}
	homog(2, 2) = 1;

	// Calculate, in both directions, the transfered points.
	NyxMat<double> hP1 = homog * points1;			  // H*P1
	NyxMat<double> invhP2 = MatInv((!homog)*homog)*(!homog)*points2;  // H\P2.

	/*
	/////
	cout << "before normalize:" << endl;
	display(points1);
	display(points2);
	display(hP1);
	display(invhP2);

	// Normalize so that the homogeneous scale 
	// parameter for all coordinates is 1.
	
	Mm newPts1 = NormalizeCrd(points1);
	Mm newPts2 = NormalizeCrd(points2);
	
	hP1 = NormalizeCrd(hP1);
	invhP2 = NormalizeCrd(invhP2);

	/////
	cout << "after normalize:" << endl;
	display(newPts1);
	display(newPts2);
	display(hP1);
	display(invhP2);
	*/

	int ptsNum = points1.GetCol();

	for (int i = 0; i < ptsNum; ++i)
	{
		double dist = 
			pow((points1(0, i) - invhP2(0, i)), 2) + 
			pow((points1(1, i) - invhP2(1, i)), 2) + 
			pow((points2(0, i) -    hP1(0, i)), 2) + 
			pow((points2(1, i) -    hP1(1, i)), 2);

		if (fabs(dist) < t)
		{
			inliers.push_back(i);
		}
	}
}


/********************************************************************
   Compute affine transform from 3 or more points.
********************************************************************/
NyxMat<double> Affine(const NyxMat<double>& points1, const NyxMat<double>& points2)
{
	int ptsNum = points1.GetCol();
	NyxMat<double> a(2*ptsNum, 6, 0);
	NyxMat<double> b(2*ptsNum, 1, 0);

	for (int i=0; i<ptsNum; ++i)
	{
		a(i*2, 0) = points1(0, i);
		a(i*2, 1) = points1(1, i);
		a(i*2, 2) = points1(2, i);

		a(i*2+1, 3) = points1(0, i);
		a(i*2+1, 4) = points1(1, i);
		a(i*2+1, 5) = points1(2, i);

		b(i*2, 0) = points2(0, i);
		b(i*2+1, 0) = points2(1, i);
	}

	// pinv(a) = (A^T * A)^(-1) * A^T
	// x=pinv(A)*b.
	NyxMat<double> x = MatInv((!a)*a)*(!a)*b;

	NyxMat<double> affine(2, 3, 0);
	// Extract affine transform.
	for (int i=0; i<2; ++i)
	{
		for (int j=0; j<3; ++j)
		{
			affine(i, j) = x(i*3+j, 0);
		}
	}

	return affine;
}

/********************************************************************
   Determine if a set of 3 pairs of matched points give rise to a 
   degeneracy in the calculation of an affine transform as needed 
   by RANSAC, by testing whether the 3 pairs are collinear.
********************************************************************/
bool IsDegenAffine(const NyxMat<double> &points1, const NyxMat<double> &points2)
{
	assert(points1.GetRow() == 3 && points2.GetRow() == 3);
	assert(points1.GetCol() == 3 && points2.GetCol() == 3);

	NyxMat<double> p11(1, 3, 0);
	p11(0, 0) = points1(0, 0);
	p11(0, 1) = points1(1, 0);
	p11(0, 2) = points1(2, 0);

	NyxMat<double> p12(1, 3, 0);
	p12(0, 0) = points1(0, 1);
	p12(0, 1) = points1(1, 1);
	p12(0, 2) = points1(2, 1);

	NyxMat<double> p13(1, 3, 0);
	p13(0, 0) = points1(0, 2);
	p13(0, 1) = points1(1, 2);
	p13(0, 2) = points1(2, 2);

	NyxMat<double> p21(1, 3, 0);
	p21(0, 0) = points2(0, 0);
	p21(0, 1) = points2(1, 0);
	p21(0, 2) = points2(2, 0);

	NyxMat<double> p22(1, 3, 0);
	p22(0, 0) = points2(0, 1);
	p22(0, 1) = points2(1, 1);
	p22(0, 2) = points2(2, 1);

	NyxMat<double> p23(1, 3, 0);
	p23(0, 0) = points2(0, 2);
	p23(0, 1) = points2(1, 2);
	p23(0, 2) = points2(2, 2);
	
	bool isDeg = 
		IsCollinear(p11, p12, p13) || 
		IsCollinear(p21, p22, p23);

	return isDeg;
}


/*------------------------------Write matrix file----------------------------*/

/********************************************************************
   Write the matrix to an output stream..
********************************************************************/
ostream & operator <<(ostream &os, const NyxMat<double> &m)
{
	int row = m.GetRow();
	int col = m.GetCol();

	for (int i=0; i<row; ++i)
	{
		for (int j=0; j<col; ++j)
		{
			os << m(i,j) << " ";
		}
		os << "\n";
	}

	return os;
}
