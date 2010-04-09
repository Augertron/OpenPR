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

/**********************************************************************
 * All rights of this (template) implemention of Class "CImgMatrix" is reserved
 * by Styx, a guy comes from Robot Vision Group in NLPR Lab, CASIA. 
 * 
 * The main goal is to provide a matlab-like interface in operating
 * image matrix. So it is a useful tool for guys be addicting in matlab but have to 
 * write c-style codes once at a time.  
 *                                                 -- Styx (11.06, 2007)
 *
 * Revised (03.04, 2009)
 *     Change class name from "CImgMatrix" to "NyxMat" for I have decided to 
 * extend this class for general use instead of only for image matrix. (And my name 
 * has been changed from Styx to Nyx too. I think Nyx is cooler 8-b ) 
 *                                                 -- Nyx (03.04, 2009)
 **************************************************************************/

#ifndef NYX_MAT_H
#define NYX_MAT_H

#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <stdexcept>

/* 圆周率 */
// const double PI = 3.1415926;

enum MatKind
{
	DIAG = 0, // diagnoal matrix
	GRIDX,    // x-directional step up matrix
	GRIDY    // y-directional step up matrix
};


template <class Type> class NyxMat
{
public:
	NyxMat();

	// gray image: depth = 1,  rgb image: depth = 3
	// construct an image matrix with default depth 1 and elements value "element"
	// only gray image class is implemented
	NyxMat(const int rows, const int cols, Type element);
		/*const unsigned int depth = 1,*/
	
	// don't initialize in construction 
	NyxMat(const int rows, const int cols); /*, const unsigned int depth = 1*/ 

	// construct 3 different kinds of matrices: diagonal, x-step-up, y-step-up
	// "nums": only for diagonal matrix, means the value of elements in the diagonal line
	NyxMat(const int rows, const int cols, MatKind whatM, Type nums = 1);
	
	// asignment between two matrices
	NyxMat& operator= (const NyxMat& sourceM);

	// copy constructor
	NyxMat(const NyxMat& sourceM);

	// crop and then copy
	NyxMat(const NyxMat& sourceM, int orix, int rows, int oriy, int cols);  
	
	//  operator overload: "()"
	Type& operator() (int i, int j);
	const Type& operator() (int i, int j) const;
	Type& operator() (int position);
	const Type& operator() (int position) const;
    
	// minus: -A
	NyxMat operator- ();                    

	// comparison between two matrics
	template <class Type2>
	NyxMat<int> operator< (const NyxMat<Type2>& compM);

	template <class Type2>
	NyxMat<int> operator<= (const NyxMat<Type2>& compM);

	// comparison between a matrix and a number
	template <class Type2>
	NyxMat<int> operator< (Type2 cmpNum);

	template <class Type2>
	NyxMat<int> operator<= (Type2 cmpNum);

	// A = A^order
	NyxMat operator^ (int order);

	// overload dot-divide and dot-multiply which runs just like as in matlab
	NyxMat DotDiv(const NyxMat& dotDivM);
	NyxMat DotMulti(const NyxMat& dotMultiM);
	NyxMat DotPower(int order);

	// extend matrix in four directions at the same time
	NyxMat ExtendImg(int hOffset, int wOffset);

	// exchange two rows (columns)
	void ExchangeRows(int rowa, int rowb);
	void ExchangeCols(int cola, int colb);

	// contract matrix in four directions at the same time (unimplemented now)
	//NyxMat ContractImg(int hOffset, int wOffset);
	
	// translation operator..........
	// by this, NyxMat<double> can easily been translated to NyxMat<int> 
	// and so does backward operation
	template <class DstType>
	operator NyxMat<DstType>() const;

	inline int GetRow() const { return m_rows; }
	inline int GetCol() const { return m_cols;  }
	inline Type** GetMatrix() const { return m_mat; }
	// inline unsigned int GetDepth() const { return depth; }

	~NyxMat();

private:
	void AllocMem(int rows, int cols);
	void DeAllocMem();

private:
	Type** m_mat;
	int m_rows;
	int m_cols;

};	


/* overload input operator "<<" and output operator ">>" */
template <class Type>
std::ostream& operator<< (std::ostream& out, const NyxMat<Type>& outImg);

template <class Type>
std::istream& operator>> (std::istream& in, NyxMat<Type>& inImg);

// overload operator "!" to implement transposition
template <class Type>
NyxMat<Type> operator! (const NyxMat<Type>& srcM);

/*
template <class Type> class Zeros : public NyxMat<Type>
{
public:
	Zeros(const unsigned int rows, const unsigned int cols);
};
 
template <class Type> class Ones : public NyxMat<Type>
{
public:
	Ones(const unsigned int rows, const unsigned int cols);
};

template <class Type> class Diag : public NyxMat<Type>
{
public:
	Diag(const unsigned int rows, const unsigned int cols, Type nums);
};
*/

// matrix addition: A + B 
template <class Type1, class Type2>
NyxMat<Type1> operator+ (const NyxMat<Type1>& srcM, const NyxMat<Type2>& addM);

// A + num
template <class Type, class TypeNum> 
NyxMat<Type> operator+ (const NyxMat<Type>& sourceM, const TypeNum nums);

// matrix substraction: A - B
template <class Type1, class Type2>
NyxMat<Type1> operator- (const NyxMat<Type1>& srcM, const NyxMat<Type2>& addM);

// A - num
template <class Type, class TypeNum> 
NyxMat<Type> operator- (const NyxMat<Type>& sourceM, const TypeNum nums);

// matrix multiplication: A * B
template <class Type1, class Type2>
NyxMat<Type1> operator* (const NyxMat<Type1>& srcM, const NyxMat<Type2>& addM);

// A * num
template <class Type, class TypeNum> 
NyxMat<Type> operator* (const NyxMat<Type>& sourceM, const TypeNum nums);

//template <class Type, class TypeNum> 
//NyxMat<Type> operator* (TypeNum nums, NyxMat<Type>& sourceM);

// matrix divition: A / B
//template <class Type1, class Type2>
//NyxMat<Type1> operator/ (const NyxMat<Type1>& srcM, const NyxMat<Type2>& addM);

template <class Type, class TypeNum> 
NyxMat<Type> operator/ (const NyxMat<Type>& sourceM, const TypeNum nums);

// exponent: only for type "double"
template <class Type>
NyxMat<double> MatExp(const NyxMat<Type>& sourceM);

// cross multiplication: 
template <class Type>
NyxMat<Type> CrossMulti(const NyxMat<Type>& vec1, const NyxMat<Type>& vec2);

// Normalization
/* which == "row",           means summerize every row, then produce a column vector,
 *                       every element of which are is sum of the row it resides
 * which == "column",        means summerize every column, then produce a row vector,
 *                       every element of which are is sum of the column it resides
 * which == "all",		     summerize all the elements in the matrix, then produce a number
 *                       of type "Type"

template <class Type>
NyxMat<Type> Sum(NyxMat<Type>& srcM, std::string which); */

template <class Type>
Type MatSum(const NyxMat<Type> &srcM);

// inverse of a square matrix
template <class Type>
NyxMat<double> MatInv(const NyxMat<Type> &A);

// get an implement of a matrix src (余子式)
template <class Type>
void GetMinor(const NyxMat<Type> &src, NyxMat<Type> &dest, int row, int col);

// calculate the determination of matrix mat
template <class Type>
Type MatDet(const NyxMat<Type> &mat);


// 2-norm of matrix
template <class Type>
double Mat2Norm(const NyxMat<Type> &srcM);

// two-dimension medieval Filters with window size of wndSize x wndSize
template <class Type>
NyxMat<Type> MedFilt2(NyxMat<Type>& srcM, int wndSize);

// a simple sortion algorithm using insertion
template <class Type>
void InsertSort(NyxMat<Type>& patch);


#include "NyxMat.cc"

#endif
