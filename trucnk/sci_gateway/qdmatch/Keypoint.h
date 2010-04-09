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

// Keypoint.h -- Declaration of the class Keypoint.

#ifndef KEYPOINT_H
#define KEYPOINT_H

/*
#include <stdlib.h>
#include <fstream.h>
#include <iostream.h>
#include <vector>

using namespace std;
*/
#include <fstream>
#include <iomanip>

#include "SiftUtil.h"


/*---------------------------------Keypoint----------------------------------*/

class Keypoint
{
public:
	Keypoint(void);
	Keypoint(double x, double y, double baseScale, double level);
	Keypoint(const Keypoint &other);
	~Keypoint(void);
	Keypoint & operator =(const Keypoint &other);

	friend bool operator ==(const Keypoint &key1, const Keypoint &key2);

	int GetDescDim(void) const;
	int GetDesc(int i) const;

	void CreateDesc(int dim, const uchar *desc);
	void DeleteDesc(void);

	//friend ostream & operator <<(ostream &os, const Keypoint &key);
	//friend istream & operator >>(istream &is, Keypoint &key);
	
public:
	double m_x; // x position of the keypoint.
	double m_y; // y position of the keypoint.
	double m_baseScale; // base scale of the keypoint.
	double m_level; // scale level of the keypoint.
	double m_scale; // scale value of the keypoint.
	double m_orien; // orientation of the keypoint.

private:
	int m_descDim; // dimension of the keypoint descriptor.
	uchar *m_desc; // local image descriptor of the keypoint.
};


int DistSquared(const Keypoint &key1, const Keypoint &key2);


ostream & operator <<(ostream &os, const Keypoint &key);
istream & operator >>(istream &is, Keypoint &key);

bool WriteKeyFile(const string& fileName, const vector<Keypoint> &keys);
bool ReadKeyFile(vector<Keypoint> &keys, const char *fileName);


#endif
