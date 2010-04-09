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

// SiftMatch.h

#ifndef SIFTMATCH_H
#define SIFTMATCH_H

#include "Keypoint.h"
#include "Match.h"
#include "KDTree.h"
#include "SiftUtil.h"


/*-------------------------------Match keypoints-----------------------------*/

// Match keypoints one by one.
int FindMatches(vector<Match> &matches, 
				const vector<Keypoint> &keyList1, 
				const vector<Keypoint> &keyList2, 
				double matchRatio);

bool CheckForMatch(int &distSq, Keypoint &matchKey, 
				   const Keypoint &key, const vector<Keypoint> &keyList, 
				   double matchRatio);


// Match by kd-tree BBF nearest neighbor searching.
int KDTreeMatch(vector<Match> &matches, 
				const vector<Keypoint> &keyList1, 
				const vector<Keypoint> &keyList2, 
				double matchRatio);


#endif
