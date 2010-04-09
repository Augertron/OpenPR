//////////////////////////////////////////////////////////////////////////////////////////////////
// Author:		Styx(Zhenhui Xu)
// Version:		0.1
// Date:		May 20, 2009 
// Description: Defination of the class Match.
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


#include "Match.h"


/*------------------------------point match class----------------------------*/

Match::Match()
{
	m_isKey = 0;
}

Match::Match(const Point &target, const Point &found, double dist) 
: m_target(target), m_found(found)
{
	m_dist = dist;
	m_isKey = 0;
}

Match::Match(const Match &other) 
: m_target(other.m_target), m_found(other.m_found)
{
	m_dist = other.m_dist;
	m_isKey = other.m_isKey;
}

//Match::~Match(void)
//{
//}

Match & Match::operator =(const Match &other)
{
	if (this == &other)
	{
		return *this;
	}

	m_target = other.m_target;
	m_found = other.m_found;
	m_dist = other.m_dist;
	m_isKey = other.m_isKey;

	return *this;
}

/*
bool Match::operator <(const Match &match) const
{
	// Compare the ZNCC correlation score.
	return (m_dist > match.m_dist);
}
*/

Point Match::GetTarget(void) const
{
	return m_target;
}

Point Match::GetFound(void) const
{
	return m_found;
}

double Match::GetDist(void) const
{
	return m_dist;
}

int Match::GetIsKey(void) const
{
	return m_isKey;
}

void Match::SetIsKey(int isKey)
{
	m_isKey = isKey;
}

/*
int Match::CompareTo(const Match &match)
{
	if (m_dist < match.m_dist)
	{
		return -1;
	}
	else if (m_dist > match.m_dist)
	{
		return 1;
	}

	return 0;
}

bool operator <(Match &match1, const Match &match2)
{
	if (match1.CompareTo(match2) == -1)
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool less_match(Match &match1, const Match &match2)
{
	if (match1.CompareTo(match2) == 1)
	{
		return true;
	}
	else
	{
		return false;
	}
}
*/


/*-------------------------Read and Write point matches----------------------*/

/********************************************************************
   Write the match to an output stream..
********************************************************************/
ostream & operator <<(ostream &os, const Match &match)
{
	//os << setiosflags(ios::fixed) << setprecision(2) 
	os	<< match.m_target.m_y << " " 
		<< match.m_target.m_x << " " 
		<< match.m_found.m_y << " " 
		<< match.m_found.m_x << " ";

	os << match.m_dist << "\n";

	return os;
}

/********************************************************************
   Read the match from an input stream.
********************************************************************/
istream & operator >>(istream &is, Match &match)
{
	is >> match.m_target.m_y >> match.m_target.m_x 
		>> match.m_found.m_y >> match.m_found.m_x 
		>> match.m_dist;

	return is;
}

/********************************************************************
   Write matches to a file.
   The file format is : 
   ----------------------------------------------
   matches number
   points1_y points1_x points2_y points2_x score
   ...
   ----------------------------------------------
********************************************************************/
bool WriteMatchFile(const string& fileName, const vector<Match> &matches)
{
	// Open the file.
	ofstream outfile(fileName.c_str());

	if (!outfile.is_open())
	{
		return false;
	}

	// Write the number of matches.
	//outfile << matches.size() << "\n";

	// Write each of the match.
	vector<Match>::const_iterator i;

	for (i=matches.begin(); i!=matches.end(); ++i)
	{
		outfile << (*i);
	}

	// Close the file.
	outfile.close();

	return true;
}


bool WriteMatchFile(const string& fileName, vector<Match>* pMatch01, const vector<Match> &match12)
{
	if (pMatch01->size() != match12.size())
	{
		cerr << "The matches in *pMatch01 and match12 should be of the same size!" << endl;
		exit(-1);
	}
	// Open the file.
	ofstream outfile(fileName.c_str());

	if (!outfile.is_open())
	{
		return false;
	}

	// Write the number of matches.
	//outfile << matches.size() << "\n";
	int i;

	for (i=0; i<match12.size(); ++i)
	{
		Match& refMatch01 = (*pMatch01)[i];
		const Match& refMatch12 = match12[i];

		outfile << refMatch01.GetTarget().GetY() << " "
			<< refMatch01.GetTarget().GetX() << " "
			<< refMatch01.GetFound().GetY() << " "
			<< refMatch01.GetFound().GetX() << " "
			<< refMatch12.GetFound().GetY() << " " 
			<< refMatch12.GetFound().GetX() << " " 
			<< endl;
	}

	// Close the file.
	outfile.close();

	return true;
}


/********************************************************************
   Read matches from a file.
********************************************************************/
bool ReadMatchFile(vector<Match> &matches, const char *fileName)
{
	// Open the file.
	ifstream file(fileName);

	if (!file.is_open())
	{
		return false;
	}

	// Read the total number of matches.
	//int n;
	//file >> n;

	// Clear the currently matches.
	matches.clear();
	
	//// Resize the vector of matches.
	//matches.resize(n);

	// Read each of the matches.
	vector<Match>::iterator i;

	Match tmp;
	while (file>>tmp)
	{
		//file >> tmp;
		matches.push_back(tmp);
		//file >> (*i++);
	}

	//for (i=matches.begin(); i!=matches.end(); ++i)
	//{
	//	file >> (*i);
	//}

	// Close the file.
	file.close();

	return true;
}


/*---------------------------------point class-------------------------------*/

Point::Point(void)
{
}

Point::Point(double x, double y)
{
	m_x = x;
	m_y = y;
}

Point::Point(const Keypoint &key)
{
	m_x = key.m_x;
	m_y = key.m_y;
}

//Point::~Point(void)
//{
//}

double Point::GetX(void) const
{
	return m_x;
}

double Point::GetY(void) const
{
	return m_y;
}
