//////////////////////////////////////////////////////////////////////////////////////////////////
// Author:		Styx(Zhenhui Xu)
// Version:		0.1
// Date:		May 20, 2009 
// Description: Heap sort algorithms.         
// Others:      Refer to the dinosaur book "Introduction to Algorithms".
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


#include "HeapSort.h"

const int INFIMIN = 0x80000000;

/*
 * Initialize a heap for sorting. For all the elements inputed have been sorted
 * increasingly, this step is relatively simple.
 * 
 * And more, handles to vecot<Match> are operated instead of themselves
 *
 */
void InitialHeap(vector<int>& heap, int size)
{
  for (int i = 0; i < size; i++)
    {
      heap.push_back(i);
    }
}


// Procedure assuring the attribution of the subheap rooted at i
void MaxHeapify(vector<Match>& match, vector<int>& heap, int i)
{
  int largest;
  int l = Left(i);
  int r = Right(i);
  int heapsize = heap.size();
  if (l < heapsize && match[heap[i]] > match[heap[l]])
    largest = l;
  else
    largest = i;

  if (r < heapsize && match[heap[largest]] > match[heap[r]])
    largest = r;

  if (largest != i)
    {
	  HeapSwap(heap, i, largest);
      MaxHeapify(match, heap, largest);
    }
}

bool BuildHeap(vector<Match>& match, vector<int>& heap)
{
  int heapsize = heap.size();
  for (int i = heapsize/2-1; i >= 0; --i)
    {
      MaxHeapify(match, heap, i);
    }
  return true;
}

/*
 * In this program, acting as a priority queue is wanted instead of sortion. 
 * So, sortion is not implemented while insertion and delete can't be omitted
 *
 */
int MaxHeapAtom(vector<int>& heap)
{
  return heap[0];
}


int HeapMaxExtract(vector<Match>& match, vector<int>& heap)
{
  int tmp;
  int heapsize = heap.size();
  if (heapsize < 1)
    {
      std::cerr << "heap is empty." << std::endl;
      exit(-1);
    }
  
  // return the max element and put the last elements in heap to head
  tmp = heap[0];
  heap[0] = heap[heapsize - 1];
  heap.erase(heap.end()-1);

  // reallign the heap and keep its attribute
  MaxHeapify(match, heap, 0);
  return tmp;
}

/*
 * Increasing key may carry controdiction to heap. So adjusting is necessary.
 * This function will be used by insertionn.
 */
void HeapIncreaseKey(vector<Match>& match, vector<int>& heap, int i, int key)
{
  if (key < heap[i])
    {
      std::cerr << "new key is smaller than current key." << std::endl;
      exit(-1);
    }

  heap[i] = key;
  int pare = Parent(i);
  while (i > 0 && match[heap[pare]] > match[heap[i]])
    {
	  HeapSwap(heap, i, pare);
      i = pare;
	  pare = Parent(pare);
    }
}

void HeapInsert(vector<Match>& match, vector<int>& heap, int key)
{
  // key will be inserted into the end of the heap firstly
  int heapsize = heap.size();
  heap.push_back(INFIMIN);
  HeapIncreaseKey(match, heap, heapsize, key);
}

void HeapSwap(vector<int>& heap, int i, int j)
{
	int tmp = heap[i];
	heap[i] = heap[j];
	heap[j] = tmp;
}
