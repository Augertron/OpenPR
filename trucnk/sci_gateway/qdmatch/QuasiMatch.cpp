//////////////////////////////////////////////////////////////////////////////////////////////////
// Author:		Styx(Zhenhui Xu)
// Version:		0.1
// Date:		May 20, 2009 
// Description: Quasi-dense propagation's wrapper.
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

#include "QuasiMatch.h"

int QuasiMatch(Image& im1, Image& im2, vector<Match>& match12, 
			   double fundF[3][3], vector<Match>* pMatch01, char* match12File)
{
	int matchNum = match12.size();
	vector<Match> corrMatches;
	corrMatches.reserve(matchNum);
	// ZNCC score threshold
	matchNum = CorrMatch(corrMatches, match12, im1, im2);

	NyxMat<double> fundM;

	int iterNum = 1;
	// If the number(matches) < REQUIRE_MATCHNUM, then it needs another propagation
	while (iterNum <= MAX_ROTATENUM && matchNum < (REQUIRE_MATCHNUM-200))
	{
		// using the propagation result of the first rotation as the new seeds
		if (iterNum != 1)
			corrMatches = match12;  

		if (pMatch01 == 0)
			matchNum = QuasiDense(fundM, match12, corrMatches, im1, im2, iterNum-1);
		else
		{
			matchNum = QuasiDense(fundM, match12, corrMatches, 
				im1, im2, iterNum-1, pMatch01, match12File);
//			printf("QuasiDense done ......\n");
			
			// in this situation, only one propagation needs, so iterNum is set to its maximum here
			iterNum = MAX_ROTATENUM;
		}

		++iterNum;
	}

	//if (iterNum == 1)    // if no propagation is applied here
	//{
	//	// fundF is no need to be changed
	//	matches = corrMatches;
	//}
	if (iterNum != 1 && matchNum > 8)
	{
		// copy the fundamental matrix
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				fundF[i][j] = fundM(i, j);
			}
		}
	}

	/////////////////////////////////// Choose the top n matches in the vector.
	//int n = 30;
	//if (n < matches.size())
	//{
	//	matches.resize(n);
	//	cout << "choose the top " << matches.size() << " matches."<< endl;
	//}

	//if (matchNum - REQUIRE_MATCHNUM >= 1500)       // too much points
	//	matches.erase(matches.begin()+ REQUIRE_MATCHNUM, matches.end());

	matchNum = match12.size();

	// WriteMatchFile(matchFile, matches);
	return matchNum;
}
