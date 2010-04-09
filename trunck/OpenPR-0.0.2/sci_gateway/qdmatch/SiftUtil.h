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

// SiftUtil.h -- common definitions and functions.

#ifndef SIFTUTIL_H
#define SIFTUTIL_H


#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cassert>
#include <ctime>
#include <vector>

#pragma warning (disable: 4786)

using namespace std;

typedef unsigned char uchar;

const double PI = 3.1415926;
const double EPS = 0.00000001;


template <class T>
inline void Swap(T &v1, T &v2)
{
	T temp = v1;
	v1 = v2;
	v2 = temp;
}

template <class T>
//inline T Max(T x, T y)
inline T Maximun(T x, T y)
{
	return (x > y) ? x : y;
}

template <class T>
//inline T Min(T x, T y)
inline T Minimun(T x, T y)
{
	return (x < y) ? x : y;
}


// Simple error handling
inline void FatalError(const char *msg)
{
	cerr << msg << endl;
	exit(1);
}


#endif
