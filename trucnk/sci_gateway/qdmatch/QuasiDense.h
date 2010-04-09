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

// QuasiDense.h

#ifndef QUASIDENSE_H
#define QUASIDENSE_H

#include <algorithm>
#include <functional>
#include <list>

#include "Image.h"
#include "Match.h"
#include "Ransac.h"
#include "Correspond.h"
#include "SiftUtil.h"
#include "HeapSort.h"
#include "NyxMat.h"

/*-----------------------------Quasi-Dense Correspondence--------------------*/

int QuasiDense(NyxMat<double> &fundM, vector<Match> &matches, 
			   const vector<Match> &corrMatches, 
			   Image &im1, Image &im2, int times,
			   vector<Match>* pMatch01 = 0, char* match12File=0);

int Propagate(vector<Match> &denseMatches, 
			  const vector<Match> &matches, NyxMat<double> &fund, double thresh,
			  Image &im1, Image &im2);

bool ** InitializeFlag(int xDim, int yDim);

//double ** InitializeGrad(Image &im);

double GradMeasure(int x, int y, Image &im);

int Resample(vector<Match> &sampledMatches, 
			 const vector<Match> &matches, 
			 Image &im1, Image &im2, int times,
			 vector<Match>* pMatch01 = 0, char* match12File=0);

// initialize confident measure for the two images
bool** InitializeConfidL(Image& img, int xDim, int yDim, double sThresh);

// find corresponding points using affine matrix, and then store that match into a vector
bool AddPnt2Vec(double x1, double y1, Image &im1, Image &im2,
				bool isKey, NyxMat<double> &affine, vector<Match>& matches);


// judge whether two points are same
bool IsSamePoint(Point& pnt1, Point& pnt2);
bool IsSamePoint(double x1, double y1, double x2, double y2);
// judge whether pnt1 is included in vector "matches"
bool IsInMatchVec(Point& pnt1, vector<Match>& matches);
// overload 
bool IsInMatchVec(int x1, int x2, vector<Match>& matches);


// assign matches into match queues of local patches
void AssignMatch2Local(vector<vector<Match> >& localMatches,
					   vector<vector<Match> >& localInterst,
					   const vector<Match>& matches, int xMargin, int yMargin, 
					   int XDim, int yDim, int xSquare, int squareNum, int square);
//overload
void AssignMatch2Local(vector<vector<Match> >& match01Pat,
					   vector<Match>* pMatch01, int xMargin, int yMargin, 
					   int XDim, int yDim, int xSquare, int squareNum, int square);

#endif
