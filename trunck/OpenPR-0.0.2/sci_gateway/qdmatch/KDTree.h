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

// KDTree.h -- Declaration of the class KDTree, HyperRect, BestEntry, HREntry; 
//          Declaration and defination of the template class SortedLimitedList.

#ifndef KDTREE_H
#define KDTREE_H


/*
#include <iostream.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <vector>
//#include <algorithm>

using namespace std;
*/
#include "Keypoint.h"
#include "SiftUtil.h"


const int VALMIN = -1;
const int VALMAX = 256;
const int INTMIN = -2147483647;
const int INTMAX = 2147483647;


template <class T>
class SortedLimitedList;

class HyperRect;
class BestEntry;
class HREntry;


/*------------------------------------KDTree-----------------------------------

  Based on "An introductory tutorial on kd-trees" by Andrew W. Moore, 
  available at http://www.ri.cmu.edu/pubs/pub_2818.html.
-----------------------------------------------------------------------------*/

class KDTree
{
public:
	KDTree(void);
	KDTree(const KDTree &other);
	~KDTree(void);
	KDTree & operator =(const KDTree &other);

	Keypoint & GetCurDom(void);
	int GetSplitDim(void) const;
	int GetDepth(void); // Compute the depth of the k-d tree.
	
	KDTree * CreateKDTree(const vector<Keypoint> &exSets);

	vector<BestEntry> NNSBBF(const Keypoint &target, int q, int searchSteps);

private:
	Keypoint ChoosePivot(int &splitDim, const vector<Keypoint> &exSets);

	Keypoint NNSBBFI(int &resDistSq, int &searchSteps, 
		int maxDistSq, int q, const Keypoint &target, 
		SortedLimitedList<BestEntry> &best, 
		SortedLimitedList<HREntry> &searchHr, HyperRect &hr);

private:
	Keypoint m_curDom; // the current element.
	int m_splitDim; // the splitting dimension for subtrees.
	KDTree *m_left; // the left k-d subtree.
	KDTree *m_right; // the right k-d subtree.
};

/*------------------------------Hyperrectangle-------------------------------*/

class HyperRect
{
	friend class KDTree;

public:
	HyperRect(void);
	HyperRect(int dim);
	HyperRect(const HyperRect &other);
	~HyperRect(void);
	HyperRect & operator =(const HyperRect &other);

	HyperRect SplitAt(int splitDim, int splitVal);
	bool IsIn(const Keypoint &target);
	bool IsInReach(const Keypoint &target, double distRad);
	double Distance(const Keypoint &target);

private:
	int m_dim; // dimension of the hyperrectangle.
	int *m_leftTop; // left top corner of the hyperrectangle.
	int *m_rightBottom; // right bottom corner of the hyperrectangle.
};

/*---------------------------------Best Entry----------------------------------

  This class is used to store the best matched keypoint and its distance 
  to the target keypoint.
-----------------------------------------------------------------------------*/

class BestEntry
{
	friend class KDTree;

public:
	BestEntry(void);
	//BestEntry(const Keypoint &neighb, double dist);
	BestEntry(const Keypoint &neighb, int distSq);//, bool squared);
	BestEntry(const BestEntry &other);
	~BestEntry(void);
	BestEntry & operator =(const BestEntry &other);

	int GetDistSq(void) const;
	//double GetDist(void) const;
	Keypoint & GetNeighb(void);

	bool Equals(const BestEntry &entry1, const BestEntry &entry2);
	int CompareTo(const BestEntry &entry);
	//friend bool operator >(BestEntry &entry1, const BestEntry &entry2);

private:
	int m_distSq; // square of the distance between matched keypoints.
	//double m_dist;
	Keypoint m_neighb; // the best matched keypoint to the target.
};

bool operator >(BestEntry &entry1, const BestEntry &entry2);

/*-------------------------Hyperrectangle Entry--------------------------------

  This class is used to store the possible searching branch in the KDTree.
-----------------------------------------------------------------------------*/

class HREntry
{
	friend class KDTree;

public:
	HREntry(void);
	HREntry(double dist, const Keypoint &pivot, 
		const HyperRect &rect, KDTree *tree);
	HREntry(const HREntry &other);
	~HREntry(void);
	HREntry & operator =(const HREntry &other);

	int CompareTo(const HREntry &entry);
	//friend bool operator >(HREntry &entry1, const HREntry &entry2);

private:
	double m_dist; // the distance between the taget keypoint to 
	               // the searching hyperrectangle.
	Keypoint m_pivot; // the pivot of 'm_tree'.
	HyperRect m_rect; // the searching hyperrectangle.
	KDTree *m_tree; // the searching subtree.
};

bool operator >(HREntry &entry1, const HREntry &entry2);


/*------------------------------SortedLimitedList------------------------------

  This class defines a sorted limited list.
  Elements in the vector are sorted increasing from left to right as : 
   m_vector[0] <= m_vector[1] <= ......
  The size of the vector is limited by 'm_maxElems'.
-----------------------------------------------------------------------------*/

template <class T>
class SortedLimitedList
{
public:
	SortedLimitedList(int maxElems);

	T & operator [](int idx);
	int GetCount(void) const; // Get the size of the vector.

	int Add(T &data); // Add a new data to the vector.
	void RemoveAt(int idx); // Remove a data in the given position.

private:
	SortedLimitedList(SortedLimitedList &other);
	SortedLimitedList & operator=(SortedLimitedList &other);

	int Set(int idx, T &data);

private:
	int m_maxElems; // max size of the vector.
	vector<T> m_vector;
};


template <class T>
SortedLimitedList<T>::SortedLimitedList(int maxElems)
{
	m_maxElems = maxElems;
}

template <class T>
T & SortedLimitedList<T>::operator [](int idx)
{
	return m_vector[idx];
}

template <class T>
int SortedLimitedList<T>::GetCount(void) const
{
	return m_vector.size();
}

template <class T>
int SortedLimitedList<T>::Add(T &data)
{
	int pos = m_vector.size();

	// Find a proper position to add 'data'.
	// Make sure the elements are in an increasing order.
	while ((pos > 0) && (m_vector[pos-1] > data))
	{
		if (pos < m_maxElems)
		{
			Set(pos, m_vector[pos-1]);
		}

		--pos;
	}

	// The size of the vector cannot exceed 'm_maxElems'.
	if (pos < m_maxElems)
	{
		Set(pos, data);
	}
	else
	{
		pos = -1;
	}

	return pos;
}

template <class T>
int SortedLimitedList<T>::Set(int idx, T &data)
{
	if (idx < m_vector.size()) // Replace a data.
	{
		m_vector[idx] = data;

		return 0;
	}
	else if (idx == m_vector.size()) // Add a new data.
	{
		m_vector.push_back(data);

		return 1;
	}
	
	return -1;
}

template <class T>
void SortedLimitedList<T>::RemoveAt(int idx)
{
	m_vector.erase(m_vector.begin() + idx);
}


#endif
