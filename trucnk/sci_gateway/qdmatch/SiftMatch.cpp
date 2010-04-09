//////////////////////////////////////////////////////////////////////////////////////////////////
// Author:		Styx(Zhenhui Xu)
// Version:		0.1
// Date:		May 20, 2009 
// Description: SIFT matching programs.
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

#include "SiftMatch.h"


/****************************************************************************
   Find matches between two keypoint lists by comparing the 
   Euclidean distance directly one by one.
   This algorithm is computational costly when the number of points is large.
*****************************************************************************/
int FindMatches(vector<Match> &matches, 
				const vector<Keypoint> &keyList1, 
				const vector<Keypoint> &keyList2, 
				double matchRatio)
{
	matches.clear();
	int count = 0;
	int distSq;
	Keypoint matchKey;

	int num = keyList1.size();
	for (int i=0; i<num; ++i)
	{
		Keypoint k(keyList1[i]);

		if (CheckForMatch(distSq, matchKey, k, keyList2, matchRatio))
		{
			//Point p1(k);
			//Point p2(matchKey);
			Match match((Point)k, (Point)matchKey, sqrt((double)distSq));
			matches.push_back(match);
			++count;
		}
	}

	return count;
}

bool CheckForMatch(int &distSq, Keypoint &matchKey, 
				   const Keypoint &key, const vector<Keypoint> &keyList, 
				   double matchRatio)
{
	Keypoint minKey;

	int dSq;
	int distSq1 = 100000000, distSq2 = 100000000;
	int num = keyList.size();

	for (int i=0; i<num; ++i)
	{
		Keypoint k(keyList[i]);

		dSq = DistSquared(key, k);

		if (dSq < distSq1)
		{
			distSq2 = distSq1;
			distSq1 = dSq;
			minKey = k;
		}
		else if (dSq < distSq2)
		{
			distSq2 = dSq;
		}
	}

	// Keep matches with (dist1 / dist2 < matchRatio).
	if (distSq1 < matchRatio * matchRatio * distSq2)
	{
		distSq = dSq;
		matchKey = minKey;
		return true;
	}
	else
	{
		return false;
	}
}

/********************************************************************
   Find matches between two keypoint lists using Best Bin First 
   k-d tree nearest neighbor searching.
********************************************************************/
int KDTreeMatch(vector<Match> &matches, 
				const vector<Keypoint> &keyList1, 
				const vector<Keypoint> &keyList2, 
				double matchRatio)
{
	matches.clear();
	int matchNum = 0;

	int maxSteps = 40;
	double matchRatioSq = matchRatio * matchRatio;

	// Build KDTree of keyList2.
#ifdef DEBUG
	cout << "build KDTree of keyList2..." << endl;
#endif

	/////
	//clock_t start = clock();

	KDTree *kd = new KDTree;
	kd = kd->CreateKDTree(keyList2);

	/////
	//clock_t finish = clock();
	//double time = (double)(finish - start) / CLOCKS_PER_SEC;
	//cout << "creating kdtree uses " << time << "s" << endl;

#ifdef DEBUG
	// Display the depth of the tree.
	int depth = kd->GetDepth();
	cout << "tree depth = " << depth << endl;

	// Display the root node.
	Keypoint root = kd->GetCurDom();
	int dim = root.GetDescDim();
	cout << "root = ";
	for (int i=0; i<dim; ++i)
	{
		cout << (int)root.GetDesc(i) << ";" ;
	}
	cout << endl;
	int splitDim = kd->GetSplitDim();
	cout << "root splitDim = " << splitDim << endl;
#endif

	// Find the nearest and second nearest neighbor 
	// for each feature in keyList1.
#ifdef DEBUG
	cout << "search in keyList2 for match of keyList1..." << endl;
#endif

	int keyNum1 = keyList1.size();

	/////
	//start = clock();

	for (int n=0; n<keyNum1; ++n)
	{
		// Input and display target.
		Keypoint target(keyList1[n]);

		// BBF search the nearest neighbors of target in the KDTree.
		vector<BestEntry> best = kd->NNSBBF(target, 2, maxSteps);

		int distSq1 = best[0].GetDistSq();
		int distSq2 = best[1].GetDistSq();
		if (distSq1 < matchRatioSq * distSq2)
		{
			Match match((Point)target, (Point)best[0].GetNeighb(), 
				sqrt((double)distSq1));
			matches.push_back(match);
			++matchNum;
		}
	}

	/////
	//finish = clock();
	//time = (double)(finish - start) / CLOCKS_PER_SEC;
	//cout << "BBFNN searching uses " << time << "s" << endl;

	// Delete KDTree.
	delete kd;

	return matchNum;
}
