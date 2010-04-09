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

#ifndef HEAPSORT_H
#define HEAPSORT_H

#include <vector>
//#include <iostream>
#include "Match.h"

//using namespace std;

/*
 * Get the index of parents and childs in left and right sub-heap.
 */
inline int Parent(int child)
{
  return (child+1)/2 - 1;
}

inline int Left(int parent)
{
  return 2*parent + 1;
}

inline int Right(int parent)
{
  return 2*parent + 2;
}


void InitialHeap(vector<int>& heap, int size);
void MaxHeapify(vector<Match>& match, vector<int>& heap, int i);

bool BuildHeap(vector<Match>& match, vector<int>& heap);
int MaxHeapAtom(vector<int>& heap);
int HeapMaxExtract(vector<Match>& match, vector<int>& heap);
void HeapInsert(vector<Match>& match, vector<int>& heap, int key);
void HeapIncreaseKey(vector<Match>& match, vector<int>& heap, int i, int key);

void HeapSwap(vector<int>& heap, int i, int j);

#endif
