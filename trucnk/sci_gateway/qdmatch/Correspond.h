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

//Correspond.h

#ifndef CORRESPOND_H
#define CORRESPOND_H

#include <algorithm>
#include <functional>

#include "Image.h"
#include "Match.h"
#include "Ransac.h"
#include "SiftUtil.h"
#include "NyxMat.h"

/*--------------------------------Correspondence-----------------------------*/

int Correspond(NyxMat<double> &F, vector<Match> &matches, 
			   const vector<Match> &initialMatches, 
			   int xDim, int yDim);


/*---------------------------------Correlation-------------------------------*/

int CorrMatch(vector<Match> &corrMatches, 
			  const vector<Match> &matches, Image &im1, Image &im2);

double NormCrossCorr(int x1, int y1, int x2, int y2, int radius, 
					 Image &im1, Image &im2);

bool ParabolaInter(double &peakPos, double &peakVal, 
				   double left, double middle, double right);


/*-----------------------Non-minimal/maximal suppression---------------------*/

int NonMinSupp(vector<Match> &reMatches, 
			   const vector<Match> &matches, 
			   int xDim, int yDim, int radius);

int NonMaxSupp(vector<Match> &reMatches, 
			   const vector<Match> &matches, 
			   int xDim, int yDim, int radius);


/*--------------------------Epipolar line constraint-------------------------*/
int EpiLineConstrn(vector<Match> &reMatches, 
				   const vector<Match> &matches, 
				   const NyxMat<double> &fund, double thresh);

int EpiLineConstrn(vector<Match> &reMatches, vector<int>& inliers,
				   const vector<Match> &matches, 
				   const NyxMat<double> &fund, double thresh);

bool IsInLimit(const Match& match, const NyxMat<double> &fund, double thresh);
bool IsInLimit(double x1, double y1,
			   double x2, double y2,
			   const NyxMat<double> &fund, double thresh);
#endif
