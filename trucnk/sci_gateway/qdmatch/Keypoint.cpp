//////////////////////////////////////////////////////////////////////////////////////////////////
// Author:		Styx(Zhenhui Xu)
// Version:		0.1
// Date:		May 20, 2009 
// Description: Definition of the class Keypoint.
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



#include "Keypoint.h"

/*-----------------------------------Keypoint--------------------------------*/

Keypoint::Keypoint(void)
{
	m_descDim = 0;
	m_desc = NULL;
}

Keypoint::Keypoint(double x, double y, double baseScale, double level)
{
	m_x = x;
	m_y = y;
	m_baseScale = baseScale;
	m_level = level;

	m_descDim = 0;
	m_desc = NULL;
}

Keypoint::Keypoint(const Keypoint &other)
{
	m_x = other.m_x;
	m_y = other.m_y;
	m_baseScale = other.m_baseScale;
	m_level = other.m_level;
	m_scale = other.m_scale;
	m_orien = other.m_orien;

	this->CreateDesc(other.m_descDim, other.m_desc);
}

Keypoint::~Keypoint(void)
{
	this->DeleteDesc();
}

Keypoint & Keypoint::operator =(const Keypoint &other)
{
	if (this == &other)
	{
		return *this;
	}

	this->DeleteDesc();
	
	m_x = other.m_x;
	m_y = other.m_y;
	m_baseScale = other.m_baseScale;
	m_level = other.m_level;
	m_scale = other.m_scale;
	m_orien = other.m_orien;

	this->CreateDesc(other.m_descDim, other.m_desc);

	return *this;
}

bool operator ==(const Keypoint &key1, const Keypoint &key2)
{
	int dim = key1.GetDescDim();
	
	for (int i=0; i<dim; ++i)
	{
		if (key1.GetDesc(i) != key2.GetDesc(i))
		{
			return false;
		}
	}

	return true;
}

int Keypoint::GetDescDim(void) const
{
	return m_descDim;
}

int Keypoint::GetDesc(int i) const
{
	assert(i>=0 && i<m_descDim);

	return (int)m_desc[i];
}

void Keypoint::CreateDesc(int dim, const uchar *desc)
{
	assert(dim >= 0);
	m_descDim = dim;

	m_desc = new uchar[m_descDim];
	if (m_desc == NULL)
	{
		FatalError("Keypoint -- Allocating memory fails!");
	}

	for (int i=0; i<m_descDim; ++i)
	{
		m_desc[i] = desc[i];
	}
}

void Keypoint::DeleteDesc(void)
{
	if (m_desc != NULL)
	{
		delete []m_desc;
	}
}


/*----------------Compute distance square between two keypoints--------------*/

int DistSquared(const Keypoint &key1, const Keypoint &key2)
{
	assert(key1.GetDescDim() == key2.GetDescDim());

	int dim = key1.GetDescDim();

	int distSq = 0;

	for (int i=0; i<dim; ++i)
	{
		int diff = key1.GetDesc(i) - key2.GetDesc(i);
		distSq += diff * diff;
	}

	return distSq;
}


/*---------------------------Read and Write features-------------------------*/

/********************************************************************
   Write the feature to an output stream.
********************************************************************/
ostream & operator <<(ostream &os, const Keypoint &key)
{
	os << setiosflags(ios::fixed) << setprecision(2) << 
		key.m_y << " " << key.m_x << " " << key.m_scale << " ";
	os << setiosflags(ios::fixed) << setprecision(3) << key.m_orien;

	int dim = key.GetDescDim();

	for (int i=0; i<dim; ++i)
	{
		// write 20 descriptor values per line
		if((i % 20) == 0)
			os << "\n";
		
		os << key.GetDesc(i) << " ";
	}
	os << "\n";

	return os;
}

/********************************************************************
   Read the feature from an input stream.
********************************************************************/
istream & operator >>(istream &is, Keypoint &key)
{
	is >> key.m_y >> key.m_x >> key.m_scale >> key.m_orien;

	int dim = 128;

	uchar *desc = new uchar[dim];
	for (int i=0; i<dim; ++i)
	{
		int tmp;
		
		is >> tmp;

		desc[i] = (uchar)tmp;
	}

	key.CreateDesc(dim, desc);

	// Modified by Styx. (Dec 13, 2007)
	// I think desc needs to be deleted.
	delete []desc;

	return is;
}

/********************************************************************
   Write keypoints to a file.
   The file format is : 
   --------------------------------------------------------
   keypoints number; descriptor dimension;
   points_y; points_x; points_scale; points_orienatition;
   points_descriptor;
   ...
   --------------------------------------------------------
********************************************************************/
bool WriteKeyFile(const string& fileName, const vector<Keypoint> &keys)
{
	// Open the file.
	ofstream outfile(fileName.c_str());

	if (!outfile.is_open())
	{
		return false;
	}

	// Write the number of keypoints and dimension of descriptor.
	outfile << keys.size() << " " << (keys[0]).GetDescDim() << endl;

	// Write each of the keypoints.
	vector<Keypoint>::const_iterator i;

	for (i=keys.begin(); i!=keys.end(); ++i)
	{
		outfile << (*i);
	}

	// Close the file.
	outfile.close();

	return true;
}

/********************************************************************
   Read keypoints from a file.
********************************************************************/
bool ReadKeyFile(vector<Keypoint> &keys, const char *fileName)
{
	// Open the file.
	ifstream file(fileName);

	if (!file.is_open())
	{
		return false;
	}

	// Read the total number of keypoints.
	int n;
	file >> n;

	// Read the dimension of each keypoint descriptor.
	// It better be 128.
	int m;
	file >> m;

	if (m != 128)
	{
		file.close();
		return false;
	}

	// Clear the currently keypoints.
	keys.clear();
	
	// Resize the vector of keypoints.
	keys.resize(n);

	// Read each of the keypoints.
	vector<Keypoint>::iterator i;

	for (i=keys.begin(); i!=keys.end(); ++i)
	{
		file >> (*i);
	}

	// Close the file.
	file.close();

	return true;
}
