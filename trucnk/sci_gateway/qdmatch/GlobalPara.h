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

#ifndef GLOBAL_PARA_H
#define GLOBAL_PARA_H

/*
 *	Some global parameters here
 */


///////////////////////////////////////
// for correlate match
///////////////////////////////////////
const int NEIGHBOR = 5;      // the range of neighbors to be searched : (-neighb~neighb)^2
const int ZNCC_WINDOW = 11;  // the rectangular correlation window size
const double SCORE_THRESH = 0.8; // the threshold to remove matches with low correlation scores
// (advice)
//        simple scene:  make SCORE_THRESH smaller. e.g, sand --> 0.7
//        complex scene: make SCORE_THRESH bigger. e.g, boat --> 0.8

const double MATCH_RATIO = 0.8;  // (nearest neighbor / secondary nearest neighbor)

const double RANSAC_THRESH = 0.001;  // thresh for ransac fit fundamental matrix

const double EPIPOLAR_CONSTRAINT_THRESH = 1.5; // thresh for eipoloar constraint

////////////////////////////////////////
// parameters for quasi_dense
////////////////////////////////////////
const double F_THRESH = 0.001;  // thresh for Ransac which is used to fit a fundamental matrix
const double E_THRESH = 1.0;    // error thresh for epipolar constraint

////////////////////////////////////////
// parameters for propagate
////////////////////////////////////////
const int PROPA_NEIGHBOR = 2;	 // radius of the neighborhood region for propagation
// (advice)
//         small image: PROPA_NEIGHBOR = 2
//         big image:   PROPA_NEIGHBOR = 3
const int PROPA_RADIUS = 2;		 // radius of the ZNCC window.
const int D_THRESH = 1;			 // disparity gradient limit.
const double C_THRESH = 0.5;	 // ZNCC value threshold.
const double S_THRESH = 0.01;	 // confidence measure threshold.


////////////////////////////////////////
//  parameters for matches' number
////////////////////////////////////////
const int REQUIRE_MATCHNUM = 1000;   // the matches needed
const int MAX_ROTATENUM = 3;           // The maximum rotation number - 1

// threshold to judge whether two points are in fact the same one
const double DIST_SAME_THRESH = 0.1;

#endif
