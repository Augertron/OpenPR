//////////////////////////////////////////////////////////////////////////////////////////////////
// Author:		Styx(Zhenhui Xu)
// Version:		0.1
// Date:		May 20, 2009 
// Description: Heap sort algorithms.         
// Others:      Defination of the class KDTree, HyperRect, BestEntry, HREntry.
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

#include "KDTree.h"


/*------------------------------------KDTree---------------------------------*/

KDTree::KDTree(void) : m_left(NULL), m_right(NULL)
{
}

KDTree::KDTree(const KDTree &other)
{
	m_curDom = other.m_curDom;
	m_splitDim = other.m_splitDim;
	m_left = NULL;
	m_right = NULL;
	if (other.m_left != NULL)
	{
		m_left = new KDTree(*other.m_left);
		if (m_left == NULL)
		{
			FatalError("KDTree -- Allocating memery fails!");
		}
	}
	if (other.m_right != NULL)
	{
		m_right = new KDTree(*other.m_right);
		if (m_right == NULL)
		{
			FatalError("KDTree -- Allocating memery fails!");
		}
	}
}

KDTree::~KDTree(void)
{
	if (m_left != NULL)
	{
		delete m_left;
		m_left = NULL;
	}
	if (m_right != NULL)
	{
		delete m_right;
		m_right = NULL;
	}
}

KDTree & KDTree::operator =(const KDTree &other)
{
	if (this == &other)
	{
		return *this;
	}

	if (m_left != NULL)
	{
		delete m_left;
		m_left = NULL;
	}
	if (m_right != NULL)
	{
		delete m_right;
		m_right = NULL;
	}

	m_curDom = other.m_curDom;
	m_splitDim = other.m_splitDim;
	m_left = NULL;
	m_right = NULL;

	if (other.m_left != NULL)
	{
		m_left = new KDTree(*other.m_left);
		if (m_left == NULL)
		{
			FatalError("KDTree -- Allocating memery fails!");
		}
	}
	if (other.m_right != NULL)
	{
		m_right = new KDTree(*other.m_right);
		if (m_right == NULL)
		{
			FatalError("KDTree -- Allocating memery fails!");
		}
	}

	return *this;
}

Keypoint & KDTree::GetCurDom(void)
{
	return m_curDom;
}

int KDTree::GetSplitDim(void) const
{
	return m_splitDim;
}

int KDTree::GetDepth(void)
{
	if (this == NULL)
	{
		return -1;
	}
	else
	{
		//return 1 + Max(m_left->GetDepth(), m_right->GetDepth());
		return 1 + Maximun(m_left->GetDepth(), m_right->GetDepth());
	}
}

/********************************************************************
   Constructing a k-d tree from a set of exemplars.
   Return the pointer of the root.
   Based on "An introductory tutorial on kd-trees" 
   by Andrew W. Moore, page 6, Table 6.3.
********************************************************************/
KDTree * KDTree::CreateKDTree(const vector<Keypoint> &exSets)
{
	int exNum = exSets.size();
	if (exNum == 0)
	{
		return NULL;
	}

	int splitDim;
	Keypoint curDom = this->ChoosePivot(splitDim, exSets);

	KDTree *cur = new KDTree;
	if (cur == NULL)
	{
		FatalError("KDTree -- Allocating memery fails!");
	}
	cur->m_curDom = curDom;
	cur->m_splitDim = splitDim;

	// Split the exemplar sets into left and right subsets 
	// relative to the splitting dimension.
	int bound = (curDom).GetDesc(splitDim);
	vector<Keypoint> leftElems, rightElems;
	
	//////////////////////////////////////////////////////
	leftElems.reserve(exNum);
	rightElems.reserve(exNum);
	//////////////////////////////////////////////////////

	for (int n=0; n<exNum; ++n)
	{
		const Keypoint &dom = exSets[n];

		// Ignore the current element.
		if (dom == curDom)
		{
			continue;
		}

		if (dom.GetDesc(splitDim) <= bound)
		{
			leftElems.push_back(dom);
		}
		else
		{
			rightElems.push_back(dom);
		}
	}

	// recursion.
	cur->m_left = this->CreateKDTree(leftElems);
	cur->m_right = this->CreateKDTree(rightElems);

	return cur;
}

/********************************************************************
   Choose a pivot and a splitting plane from which to build 
   the root of a k-d tree.
   Here we choose a pivot which splits the exemplar set in the 
   middle of the range of the most spread dimension.
********************************************************************/
Keypoint KDTree::ChoosePivot(int &splitDim, const vector<Keypoint> &exSets)
{
	int exNum = exSets.size();
	int dim = exSets[0].GetDescDim();
	int i, n;

	// Find each dimension's minimum and maximum elements for every exemplar.
	int *dimMin = new int[dim];
	if (dimMin == NULL)
	{
		FatalError("KDTree::ChoosePivot() -- Allocating memory fails!");
	}
	int *dimMax = new int[dim];
	if (dimMax == NULL)
	{
		FatalError("KDTree::ChoosePivot() -- Allocating memory fails!");
	}
	for (i=0; i<dim; ++i)
	{
		dimMin[i] = VALMAX;
		dimMax[i] = VALMIN;
	}

	for (n=0; n<exNum; ++n)
	{
		const Keypoint &dom = exSets[n];

		for (i=0; i<dim; ++i)
		{
			int dimElem = dom.GetDesc(i);

			if (dimElem < dimMin[i])
			{
				dimMin[i] = dimElem;
			}
			if (dimElem > dimMax[i])
			{
				dimMax[i] = dimElem;
			}
		}
	}

	// Find the maximum range dimension.
	int *dimDiff = new int[dim];
	if (dimDiff == NULL)
	{
		FatalError("KDTree::ChoosePivot() -- Allocating memory fails!");
	}
	int maxDiff = 0;
	int maxDiffDim = 0;

	for (i=0; i<dim; ++i)
	{
		dimDiff[i] = dimMax[i] - dimMin[i];

		if (dimDiff[i] > maxDiff)
		{
			maxDiff = dimDiff[i];
			maxDiffDim = i;
		}
	}

	// The splitting dimension is maxDiffDim, 
	// now find a exemplar as close to the arithmetic middle as possible.
	double middle = (maxDiff / 2.0) + dimMin[maxDiffDim];
	double exMinDiff = VALMAX;
	Keypoint exemplar;

	for (n=0; n<exNum; ++n)
	{
		const Keypoint &dom = exSets[n];
		double curDiff = fabs(dom.GetDesc(maxDiffDim) - middle);

		if (curDiff < exMinDiff)
		{
			exMinDiff = curDiff;
			exemplar = dom;
		}
	}

	/*
	// The splitting dimension is maxDiffDim,
	// now find the median exemplar to this dimension.
	vector<int> dimValue;

	for (n=0; n<exNum; ++n)
	{
		dimValue.push_back(exSets[n].GetDesc(maxDiffDim));
	}

	sort(dimValue.begin(), dimValue.end());
	int median = dimValue[exNum/2];
	Keypoint exemplar;

	for (n=0; n<exNum; ++n)
	{
		const Keypoint &dom = exSets[n];
		if (dom.GetDesc(maxDiffDim) == median)
		{
			exemplar = dom;
			break;
		}
	}
	*/

	delete []dimMin;
	delete []dimMax;
	delete []dimDiff;

	// Return the values.
	splitDim = maxDiffDim;

	return exemplar;
}

/********************************************************************
   Best Bin First Nearest Neighbor Searching.
   Return the best matched keypoints to the given target.
   'q' -- we found 'q' best matches for the given target; 
          'q' = 2 means the nearest neighbor and the second nearest 
		  neighbor are found.
   'searchSteps' -- limit the number of leaf nodes to be examined.

   The Nearest Neighbor Algorithm is based on 
   "An introductory tutorial on kd-trees" by Andrew W. Moore, 
   page 8, Table 6.4.
   The Best Bin First Algorithm is based on 
   "Shape indexing using approximate nearest-neighbour search in 
   high-dimensional spaces" by Jeffrey S. Beis and David G. Lowe.
********************************************************************/
vector<BestEntry> KDTree::NNSBBF(const Keypoint &target, 
								 int q, int searchSteps)
{
	HyperRect hr(target.GetDescDim());

	SortedLimitedList<BestEntry> best(q);
	SortedLimitedList<HREntry> searchHr(searchSteps);

	/////////////////////////////////////////////////////////////
	//clock_t start = clock();
	/////////////////////////////////////////////////////////////

	int dummyDist;
	Keypoint nearest = NNSBBFI(dummyDist, searchSteps, 
		INTMAX, q, target, best, searchHr, hr);

	/////////////////////////////////////////////////////////////
	//clock_t finish = clock();
	//double time = (double)(finish - start) / CLOCKS_PER_SEC;
	//cout << "time = " << time << endl;
	/////////////////////////////////////////////////////////////

	vector<BestEntry> bestEntry(q);
	for (int i=0; i<q; ++i)
	{
		bestEntry[i] = best[i];
	}

	return bestEntry;
}

Keypoint KDTree::NNSBBFI(int &resDistSq, int &searchSteps, 
						 int maxDistSq, int q, const Keypoint &target, 
						 SortedLimitedList<BestEntry> &best, 
						 SortedLimitedList<HREntry> &searchHr, 
						 HyperRect &hr)
{
	resDistSq = INTMAX;
	Keypoint pivot = m_curDom;
	int splitDim = m_splitDim;
	int splitVal = pivot.GetDesc(splitDim);

	BestEntry be(pivot, DistSquared(target, pivot));
	best.Add(be);

	/////////////////////////////////////////////////////////////
	//clock_t start = clock();
	/////////////////////////////////////////////////////////////

	HyperRect leftHr = hr;
	HyperRect rightHr = leftHr.SplitAt(splitDim, splitVal);

	HyperRect nearerHr;
	HyperRect furtherHr;
	KDTree *nearerKd = NULL;
	KDTree *furtherKd = NULL;
	
	// Determine which child contains the 'target'.
	if (target.GetDesc(splitDim) <= splitVal)
	{
		nearerHr = leftHr;
		furtherHr = rightHr;
		nearerKd = m_left;
		furtherKd = m_right;
	}
	else
	{
		nearerHr = rightHr;
		furtherHr = leftHr;
		nearerKd = m_right;
		furtherKd = m_left;
	}

	/////////////////////////////////////////////////////////////
	//clock_t finish = clock();
	//double time = (double)(finish - start) / CLOCKS_PER_SEC;
	//cout << "time=" << time << ";";
	/////////////////////////////////////////////////////////////

	// Search for the initial child.
	Keypoint nearest;
	int distSq;

	HREntry he(furtherHr.Distance(target), pivot, furtherHr, furtherKd);
	searchHr.Add(he);

	// Bottom reached.
	if (nearerKd == NULL)
	{
		distSq = INTMAX;
	}
	else
	{
		nearest = nearerKd->NNSBBFI(distSq, searchSteps, 
			maxDistSq, q, target, best, searchHr, nearerHr);
	}

	searchSteps -= 1;
	if (searchSteps <= 0)
	{
		resDistSq = distSq;
		
		return nearest;
	}

	// Restrict the maximum radius in which 
	// any possible closer point could lie.
	if (best.GetCount() >= q)
	{
		maxDistSq = best[q-1].m_distSq;
	}
	else
	{
		maxDistSq = INTMAX;
	}

	HREntry hre;
	if (searchHr.GetCount() > 0)
	{
		hre = searchHr[0];
		searchHr.RemoveAt(0);

		furtherHr = hre.m_rect;
		furtherKd = hre.m_tree;
		pivot = hre.m_pivot;
	}

	// Check whether there is any space in the hyperrectangle 
	// of the further child which lies within this radius.
	if (furtherHr.IsInReach(target, sqrt(static_cast<double>(maxDistSq))))
	{
		int ptDistSq = DistSquared(pivot, target);
		if (ptDistSq < distSq)
		{
			nearest = pivot;
			distSq = ptDistSq;
			maxDistSq = distSq;
		}

		int tempDistSq;
		Keypoint tempNearest;
		if (furtherKd == NULL)
		{
			tempDistSq = INTMAX;
		}
		else
		{
			tempNearest = furtherKd->NNSBBFI(tempDistSq, searchSteps, 
				maxDistSq, q, target, best, searchHr, furtherHr);
		}

		if (tempDistSq < distSq)
		{
			nearest = tempNearest;
			distSq = tempDistSq;
		}
	}

	resDistSq = distSq;

	return nearest;
}

/*----------------------------------HyperRect--------------------------------*/

HyperRect::HyperRect(void)
{
	m_leftTop = NULL;
	m_rightBottom = NULL;
}

HyperRect::HyperRect(int dim)
{
	m_dim = dim;

	m_leftTop = new int[m_dim];
	if (m_leftTop == NULL)
	{
		FatalError("HyperRect -- Allocating memery fails!");
	}

	m_rightBottom = new int[m_dim];
	if (m_rightBottom == NULL)
	{
		FatalError("HyperRect -- Allocating memery fails!");
	}

	for (int i=0; i<m_dim; ++i)
	{
		m_leftTop[i] = INTMIN;
		m_rightBottom[i] = INTMAX;
	}
}

HyperRect::HyperRect(const HyperRect &other)
{
	m_dim = other.m_dim;

	m_leftTop = new int[m_dim];
	if (m_leftTop == NULL)
	{
		FatalError("HyperRect -- Allocating memery fails!");
	}

	m_rightBottom = new int[m_dim];
	if (m_rightBottom == NULL)
	{
		FatalError("HyperRect -- Allocating memery fails!");
	}

	for (int i=0; i<m_dim; ++i)
	{
		m_leftTop[i]= other.m_leftTop[i];
		m_rightBottom[i] = other.m_rightBottom[i];
	}
}

HyperRect::~HyperRect(void)
{
	if (m_leftTop != NULL)
	{
		delete []m_leftTop;
	}

	if (m_rightBottom != NULL)
	{
		delete []m_rightBottom;
	}
}

HyperRect & HyperRect::operator =(const HyperRect &other)
{
	if (this == &other)
	{
		return *this;
	}

	if (m_leftTop != NULL)
	{
		delete []m_leftTop;
	}

	if (m_rightBottom != NULL)
	{
		delete []m_rightBottom;
	}

	m_dim = other.m_dim;

	m_leftTop = new int[m_dim];
	if (m_leftTop == NULL)
	{
		FatalError("HyperRect -- Allocating memery fails!");
	}

	m_rightBottom = new int[m_dim];
	if (m_rightBottom == NULL)
	{
		FatalError("HyperRect -- Allocating memery fails!");
	}

	for (int i=0; i<m_dim; ++i)
	{
		m_leftTop[i]= other.m_leftTop[i];
		m_rightBottom[i] = other.m_rightBottom[i];
	}

	return *this;
}

/********************************************************************
   Split the hyperrectangle at the dimension 'splitDim' 
   by the value 'splitVal'.
********************************************************************/
HyperRect HyperRect::SplitAt(int splitDim, int splitVal)
{
	assert(m_leftTop[splitDim] < splitVal || 
		m_rightBottom[splitDim] >= splitVal);

	HyperRect rect(*this);
	m_rightBottom[splitDim] = splitVal;
	rect.m_leftTop[splitDim] = splitVal;

	return rect;
}

/********************************************************************
   Check whether the point 'target' is in the hyperrectangle.
********************************************************************/
bool HyperRect::IsIn(const Keypoint &target)
{
	assert(target.GetDescDim() == m_dim);

	for (int i=0; i<m_dim; ++i)
	{
		int elem = target.GetDesc(i);
		if (elem < m_leftTop[i] || elem >= m_rightBottom[i])
		{
			return false;
		}
	}

	return true;
}

/********************************************************************
   Check to see if the hyperrectangle intersects with a hypersphere 
   radius 'distRad' centered at point 'target'.
********************************************************************/
bool HyperRect::IsInReach(const Keypoint &target, double distRad)
{
	return (Distance(target) < distRad);
}

/********************************************************************
   Find the point in the hyperrectangle which is closest to 
   'target', the distance is between them.
********************************************************************/
double HyperRect::Distance(const Keypoint &target)
{
	int closestPointN;
	int distance = 0;

	for (int i=0; i<m_dim; ++i)
	{
		int elem = target.GetDesc(i);
		int hrMin = m_leftTop[i];
		int hrMax = m_rightBottom[i];

		if (elem <= hrMin)
		{
			closestPointN = hrMin;
		}
		else if (elem > hrMin && elem < hrMax)
		{
			closestPointN = elem;
		}
		else //if (elem >= hrMax)
		{
			closestPointN = hrMax;
		}

		int dimDist = elem - closestPointN;
		distance += dimDist * dimDist;
	}

	return sqrt((double)distance);
}

/*----------------------------------BestEntry--------------------------------*/

BestEntry::BestEntry(void)
{
}

/*
BestEntry::BestEntry(const Keypoint &neighb, double dist) 
: m_neighb(neighb)
{
	m_dist = dist;
}
*/

BestEntry::BestEntry(const Keypoint &neighb, int distSq)//, bool squared) 
: m_neighb(neighb)
{
	m_distSq = distSq;
}

BestEntry::BestEntry(const BestEntry &other) 
: m_neighb(other.m_neighb)
{
	//m_dist = other.m_dist;
	m_distSq = other.m_distSq;
}

BestEntry::~BestEntry(void)
{
}

BestEntry & BestEntry::operator =(const BestEntry &other)
{
	if (this == &other)
	{
		return *this;
	}

	//m_dist = other.m_dist;
	m_distSq = other.m_distSq;
	m_neighb = other.m_neighb;

	return *this;
}

int BestEntry::GetDistSq(void) const
{
	return m_distSq;
}

/*
double BestEntry::GetDist(void) const
{
	return m_dist;
}
*/

Keypoint & BestEntry::GetNeighb(void)
{
	return m_neighb;
}

bool BestEntry::Equals(const BestEntry &entry1, const BestEntry &entry2)
{
	return (entry1.m_neighb == entry2.m_neighb);
}

int BestEntry::CompareTo(const BestEntry &entry)
{
	if (m_distSq < entry.m_distSq)
	{
		return -1;
	}
	else if (m_distSq > entry.m_distSq)
	{
		return 1;
	}

	return 0;
}

bool operator >(BestEntry &entry1, const BestEntry &entry2)
{
	if (entry1.CompareTo(entry2) == 1)
	{
		return true;
	}
	else
	{
		return false;
	}
}

/*------------------------------------HREntry--------------------------------*/

HREntry::HREntry(void) : m_tree(NULL)
{
}

HREntry::HREntry(double dist, const Keypoint &pivot, 
				 const HyperRect &rect, KDTree *tree) 
				 : m_pivot(pivot), m_rect(rect)
{
	m_dist = dist;
	m_tree = tree;
}

HREntry::HREntry(const HREntry &other) 
: m_pivot(other.m_pivot), m_rect(other.m_rect)
{
	m_dist = other.m_dist;
	m_tree = other.m_tree;
}

HREntry::~HREntry(void)
{
}

HREntry & HREntry::operator =(const HREntry &other)
{
	if (this == &other)
	{
		return *this;
	}

	m_dist = other.m_dist;
	m_pivot = other.m_pivot;
	m_rect = other.m_rect;
	m_tree = other.m_tree;

	return *this;
}

int HREntry::CompareTo(const HREntry &entry)
{
	if (m_dist < entry.m_dist)
	{
		return -1;
	}
	else if (m_dist > entry.m_dist)
	{
		return 1;
	}

	return 0;
}

bool operator >(HREntry &entry1, const HREntry &entry2)
{
	if (entry1.CompareTo(entry2) == 1)
	{
		return true;
	}
	else
	{
		return false;
	}
}
