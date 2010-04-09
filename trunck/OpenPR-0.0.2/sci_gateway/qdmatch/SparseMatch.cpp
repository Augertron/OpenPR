//////////////////////////////////////////////////////////////////////////////////////////////////
// Author:		Styx(Zhenhui Xu)
// Version:		0.1
// Date:		May 20, 2009 
// Description: Sparse matching's wrapper.
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

#include "SparseMatch.h"

/*
 * Extracting feature points from two corresponding images and then matching them
 * Arguments:
 *		im1, im2 -- two images represented by class Image
 *      matches  -- match points computed
 *
 * Copyright(c)Xu Zhenhui, Robot Vision Group, NLPR
 * March 2nd, 2008
 * Modified in July 23th, 2008
 *
 */
int SparseMatch(Image &im1, Image &im2, vector<Match>& matches, double fundF[3][3])
{
	//extract SIFT feature in two input images.
	double constrastThresh = 0.01;

	// Extract SIFT features for image1.
	vector<Keypoint> keyList1 = SIFT(im1, constrastThresh);

	// Extract SIFT features for image2.
	vector<Keypoint> keyList2 = SIFT(im2, constrastThresh);

	// // Read key file1
	//cout << "input keyList1 ..." << endl;
	//vector<Keypoint> keyList1;
	//ReadKeyFile(keyList1, keyFile1);
	//cout << keyList1.size() << " keys in image1." << endl;
	//
	//// Read key file2.
	//cout << "input keyList2 ..." << endl;
	//vector<Keypoint> keyList2;
	//ReadKeyFile(keyList2, keyFile2);
	//cout << keyList2.size() << " keys in image2." << endl;
	
	// Compute correspondences between two images.
	int matchNum = 0;
	
	// Get initial matches from extracted features.
	double matchRatio = MATCH_RATIO;
	vector<Match> initialMatches; 
	//initialMatches.reserve(Min(keyList1.size(), keyList2.size())/2);
	initialMatches.reserve(Minimun(keyList1.size(), keyList2.size())/2);

	matchNum = KDTreeMatch(initialMatches, keyList1, keyList2, matchRatio);

	vector<Match> firstMatches;
	int radius = 5;

	// non-minimum (distance) suppression
	matchNum = NonMinSupp(firstMatches, initialMatches, im1.GetXDim(), im2.GetYDim(), radius);

	//vector<Match> matches;
	matches.reserve(matchNum);
	NyxMat<double> F;  // fundamental matrix needs to put out

	// epipolar constraint
	matchNum = Correspond(F, matches, firstMatches, im1.GetXDim(), im2.GetYDim());

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			fundF[i][j] = F(i, j);
		}
	}

	return matchNum;
	// This part will be moved into propagation model
	//// Compute correlation scores for current matches.
	//vector<Match> corrMatches;
	//corrMatches.reserve(matchNum);

	//matchNum = CorrMatch(corrMatches, matches, im1, im2);
	//cout << "---------------------------" << 
	//	matchNum << " matches after correlation." << endl;

	////WriteMatchFile(matchFile, corrMatches);
}
