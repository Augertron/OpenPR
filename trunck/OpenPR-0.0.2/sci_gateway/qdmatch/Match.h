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

// Match.h -- Declaration of the class Match.

#ifndef MATCH_H
#define MATCH_H


#include "Keypoint.h"


class Match;

ostream & operator <<(ostream &os, const Match &match);
istream & operator >>(istream &is, Match &match);


class Point
{
public:
	Point(void);
	Point(double x, double y);
	Point(const Keypoint &key);
	// ~Point(void);

	double GetX(void) const;
	double GetY(void) const;

	friend ostream & operator <<(ostream &os, const Match &match);
	friend istream & operator >>(istream &is, Match &match);

private:
	double m_x;
	double m_y;
};


class Match
{
public:
	Match();
	Match(const Point &target, const Point &found, double dist);
	Match(const Match &other);
	//~Match(void);
	Match & operator =(const Match &other);

	friend ostream & operator <<(ostream &os, const Match &match);
	friend istream & operator >>(istream &is, Match &match);

	Point GetTarget(void) const;
	Point GetFound(void) const;
	double GetDist(void) const;
	
	int GetIsKey(void) const;
	void SetIsKey(int isKey);

	inline bool operator< (const Match& match) const
	{
		// Compare the ZNCC correlation score.
		return (m_dist > match.m_dist);
	}

	inline bool operator> (const Match& match) const
	{
		// Compare the ZNCC correlation score.
		return (m_dist <= match.m_dist);
	}
	//int CompareTo(const Match &match);

private:
	Point m_target; // matched point from the first image.
	Point m_found; // matched point from the second image.

	// Before correlation, it is distance between the matched points.
	// After correlation, it is correlation score between the match.
	double m_dist;

	// Mark if the match includes an interest point, 1 yes, 0 no.
	int m_isKey;
};


//bool operator <(Match &match1, const Match &match2);

// Sometimes m_dist stores the similarity measure value like 
// correlation coefficient, which is the larger, the better.
//bool less_match(Match &match1, const Match &match2);


//ostream & operator <<(ostream &os, const Match &match);
//istream & operator >>(istream &is, Match &match);

bool WriteMatchFile(const string& fileName, const vector<Match> &matches);
bool ReadMatchFile(vector<Match> &matches, const char *fileName);
bool WriteMatchFile(const string& fileName, vector<Match>* pMatch01, const vector<Match> &match12);
#endif


