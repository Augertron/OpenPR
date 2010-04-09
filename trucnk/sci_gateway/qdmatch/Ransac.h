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

#ifndef RANSAC_H
#define RANSAC_H

#include "NyxMat.h"
#include "Match.h"
#include "SiftUtil.h"


//const double EPS = 0.00000001;

ostream & operator <<(ostream &os, const NyxMat<double> &m);

enum FitModel
{
    eFundmatrix = 0,  // fundamental matrix model.
    eHomography = 1,  // homography model.
	eAffine     = 2   // affine model.
};


/*------------------------RANSAC fitting general model-----------------------*/

int Ransac(NyxMat<double> &m, vector<int> &inliers, 
		   const NyxMat<double> &points1, const NyxMat<double> &points2, 
		   int n, double t, FitModel fit, int maxTrials=1000, int maxDataTrials=100);

void DistFun(vector<int> &inliers, const NyxMat<double> &m, 
			 const NyxMat<double> &points1, const NyxMat<double> &points2, 
			 double t, FitModel fit);

NyxMat<double> FittingFun(const NyxMat<double> &points1, 
						  const NyxMat<double> &points2, FitModel fit);

bool DegenFun(const NyxMat<double> &points1, const NyxMat<double> &points2, FitModel fit);

/*------------------RANSAC fitting fundamental matrix model------------------*/

int RansacFitFund(NyxMat<double> &model, vector<int> &inliers, 
				  const vector<Match> &matches, double t);

void FundDist(vector<int> &inliers, const NyxMat<double> &fundM, 
			  const NyxMat<double> &points1, const NyxMat<double> &points2, double t);

NyxMat<double> Fundmatrix(const NyxMat<double> &points1, const NyxMat<double> &points2);

/*--------------------------------Normalization------------------------------*/

NyxMat<double> NormalizePts(NyxMat<double> &trans, const NyxMat<double> &pts);

NyxMat<double> NormalizeCrd(const NyxMat<double> &pts);


/*-----------------------RANSAC fitting homography model---------------------*/
int RansacFitHomog(NyxMat<double> &model, vector<int> &inliers, 
				   const vector<Match> &matches, double t);

void HomogDist(vector<int> &inliers, const NyxMat<double> &homog, 
			   const NyxMat<double> &points1, const NyxMat<double> &points2, double t);

NyxMat<double> Homography(const NyxMat<double>& points1, const NyxMat<double>& points2);

bool IsDegenHomog(const NyxMat<double>& points1, const NyxMat<double> &points2);

bool IsCollinear(const NyxMat<double> &p1, 
				 const NyxMat<double> &p2, 
				 const NyxMat<double> &p3);

/*-----------------------RANSAC fitting affine model---------------------*/
int RansacFitAffine(NyxMat<double> &model, vector<int> &inliers, 
					const vector<Match> &matches, double t, 
					int maxTrials=1000, int maxDataTrials=100);

void AffineDist(vector<int> &inliers, const NyxMat<double> &affine, 
				const NyxMat<double> &points1, const NyxMat<double> &points2, double t);

NyxMat<double> Affine(const NyxMat<double>& points1, const NyxMat<double>& points2);

bool IsDegenAffine(const NyxMat<double> &points1, const NyxMat<double> &points2);

#endif
