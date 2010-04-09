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

/*************************************************************************
*                  implementation of Class "NyxMat"
*************************************************************************/

#ifndef NYX_MAT_CC
#define NYX_MAT_CC

// #include "NyxMat.h"
#include <iostream>
#include <cstdlib>  // for exit(...)

template <class Type>
void NyxMat<Type>::AllocMem(int rows, int cols)
{
	try	{
		if (cols <= 0 && rows <= 0)
			throw std::out_of_range("AllocMem error in NyxMat: rows or cols must be positive.");

		m_mat = new Type*[rows];
		if (m_mat == 0)
			throw std::bad_alloc();
		
		for (int i = 0; i < rows; i++)
		{
			m_mat[i] = new Type[cols];
			if (m_mat[i] == 0)
				throw std::bad_alloc();
		}

	} 
	catch (const std::out_of_range& e) 
	{
		std::cerr << e.what() << std::endl;
		throw;   // re-throw this exception to its provoker
	} 
	catch (const std::bad_alloc& e) 
	{
		this->DeAllocMem();
		std::cerr << "NyxMat: bad_alloc in AllocMem." << std::endl;
		throw;  // re-throw this exception to its provoker
	}
}


template <class Type>
void NyxMat<Type>::DeAllocMem()
{
	if (m_mat != 0)
	{
		for (int i = 0; i < m_rows; i++)
		{
			delete []m_mat[i];
		}
		delete []m_mat;
	}
}



template <class Type>
NyxMat<Type>::NyxMat()
{
	m_rows = 0;
	m_cols = 0;
	m_mat = 0;
}


template <class Type>
NyxMat<Type>::NyxMat(const int rows, const int cols, Type element)
	: m_rows(rows), m_cols(cols)
{
	this->AllocMem(rows, cols);
	// initialization
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			m_mat[i][j] = element;
}


template <class Type>
NyxMat<Type>::NyxMat(const int rows, const int cols)
	: m_rows(rows), m_cols(cols)
{
	// no initialization, only allocate memory
	this->AllocMem(rows, cols);
}



template <class Type>
NyxMat<Type>::NyxMat(const int rows, const int cols, MatKind whatM, Type nums)
	: m_rows(rows), m_cols(cols)
{
	// rows should equal with cols for it is a diagonal matrix
	if ( (whatM == DIAG) && rows != cols)	
	{
		std::cout << "NyxMat<Type>::NyxMat(rows, cols, DIAG, nums)'s WARNING:" << std::endl;
		std::cout << "	rows != cols in making a diagonal matrix." << std::endl;
		std::cout << "	It's queer, so make sure it is what you want." << std::endl;	
	}

	
	switch(whatM)
	{
	case(DIAG):
		// for diagonal matrix
		this->AllocMem(rows, cols);
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
			{
				if (i == j)
					m_mat[i][j] = nums;
				else
					m_mat[i][j] = 0;
			}
		}
		break;

	case(GRIDX):
		// produce x-direction step-up matrix like matlab's meshgrid does
		this->AllocMem(rows, cols);
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
			{
				m_mat[i][j] = i;
			}
		}
		break;

	case(GRIDY):   
		// produce y-direction step-up matrix like matlab's meshgrid does
		this->AllocMem(rows, cols);
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
			{
				m_mat[i][j] = j;
			}
		}
		break;
	default:
		std::cerr << "NyxMat: invalid argument,"
			<< " please use DIAG, GRIDX or GRIDY instead" << std::endl;
		exit(2);
		break;
	}
}



template <class Type>
NyxMat<Type>& NyxMat<Type>::operator= (const NyxMat<Type>& srcM)
{
	//if ( (cols != srcM.GetCol()) || (rows != srcM.GetRow()) )
	//	throw SizeMismatch();
	
	if (this != &srcM)           // make sure it is not a self-assignment
	{
		this->DeAllocMem();
		this->m_rows = srcM.GetRow();
		this->m_cols = srcM.GetCol();
		this->AllocMem(this->m_rows, this->m_cols);

		for (int i = 0; i < this->m_rows; i++)
			for (int j = 0; j < this->m_cols; j++)
				m_mat[i][j] = srcM(i, j);
	}

	return *this;
}


// copy constructor
template <class Type>
NyxMat<Type>::NyxMat(const NyxMat<Type>& srcM)
	: m_rows(srcM.GetRow()), m_cols(srcM.GetCol()) 
{
	this->AllocMem(m_rows, m_cols);

	for (int i = 0; i < m_rows; i++)
		for (int j = 0; j < m_cols; j++)
			m_mat[i][j] = srcM(i, j);
}


// crop and copy -- construct a matrix using the region which starts from point(oriR, oriC),
// with size (offsetR x offsetC) in a matrix "srcM"
template <class Type>
NyxMat<Type>::NyxMat(const NyxMat<Type>& srcM, int oriR, int offsetR, 
						   int oriC, int offsetC)
{
	if ( (oriR < 0) || (oriC < 0) || 
		((oriR+offsetR) > srcM.GetRow()) || ((oriC+offsetC) > srcM.GetCol()))
	{
		std::cerr << "NyxMat: Sub region must be inside of the matrix." << std::endl;
		exit(2);
	}
	
	m_rows = offsetR;
	m_cols = offsetC;
	this->AllocMem(m_rows, m_cols);

	for (int i = 0; i < m_rows; i++)
	{
		for (int j = 0; j < m_cols; j++)
		{
			m_mat[i][j] = srcM(oriR+i, oriC+j);
		}
	}
}


//  overload operator "()", then the class operates like a matrix in matlab
template <class Type>
Type& NyxMat<Type>::operator() (int i, int j)
{
	// i -- rows, j -- cols
	if ( (i >= m_rows) || (j >= m_cols) || (i < 0) || (j < 0) )
	{
		std::cerr << "Nyxmat::operator(i, j): out of range." << std::endl;
		exit(2);
	}

	return m_mat[i][j];
}


template <class Type>
const Type& NyxMat<Type>::operator() (int i, int j) const
{
	// i -- rows, j -- cols
	if ( (i >= m_rows) || (j >= m_cols) || (i < 0) || (j < 0) )
	{
		std::cerr << "Nyxmat::operator(i, j): out of range." << std::endl;
		exit(2);
	}

	return m_mat[i][j];
}


// In fact, these two functions are useless
template <class Type>
Type& NyxMat<Type>::operator() (int position)
{
	if ( position >= m_rows*m_cols || position < 0)
	{
		std::cerr << "Nyxmat::operator(i): out of range." << std::endl;
		exit(2);
	}

	return m_mat[position/m_cols][position%m_cols];
}


template <class Type>
const Type& NyxMat<Type>::operator() (int position) const
{
	if ( position >= m_rows*m_cols || position < 0)
	{
		std::cerr << "Nyxmat::operator(i): out of range." << std::endl;
		exit(2);
	}

	return m_mat[position/m_cols][position%m_cols];
}


// exchange two rows of matrix
template <class Type>
void NyxMat<Type>::ExchangeRows(int rowa, int rowb)
{
	Type tmp;

	if (rowa<0 || rowa>=m_rows || rowb<0 || rowb>=m_rows)
	{
		std::cout << "ExchangeRows ERROR: rowa or rowb must be belong to [0, m_rows)" << std::endl;
		exit(1);
	}

	for (int i = 0; i < m_cols; ++i)
	{
		tmp = m_mat[rowa][i];
		m_mat[rowa][i] = m_mat[rowb][i];
		m_mat[rowb][i] = tmp;
	}
}

// exchange two columns of matrix
template <class Type>
void NyxMat<Type>::ExchangeCols(int cola, int colb)
{
	Type tmp;

	if (cola<0 || cola>=m_cols || colb<0 || colb>=m_cols)
	{
		std::cout << "ExchangeCols ERROR: cola or colb must be belong to [0, m_cols)" << std::endl;
		exit(1);
	}

	for (int i = 0; i < m_rows; ++i)
	{
		tmp = m_mat[i][cola];
		m_mat[i][cola] = m_mat[i][colb];
		m_mat[i][colb] = tmp;
	}
}

// overload operator "!" to implement transposition
template <class Type>
NyxMat<Type> operator! (const NyxMat<Type>& srcM)
{
	int rows = srcM.GetRow();
	int cols = srcM.GetCol();
	NyxMat<Type> tmpM(cols, rows);             // A' is n x m if A is m x n
	for ( int i = 0; i < cols; i++)
	{
		for ( int j = 0; j < rows; j++)
		{
			tmpM(i, j) = srcM(j, i);
		}
	}
	return tmpM;
}


// A^2 = A * A
template <class Type>
NyxMat<Type> NyxMat<Type>::operator^ (int order)
{
	if (order < 0)
	{
		std::cerr << " Temporarily, only positive integer orders are supported!" << std::endl;
		std::cerr << " Program will quit........." << std::endl;
		exit(3);
	}
	

	if (order == 0)
		return NyxMat<Type>(m_rows, m_cols, 1);

	// a diagonal matrix with elements 1
	NyxMat<Type> tmpM(m_rows, m_cols, DIAG, 1);

	for (int i = 0; i < order; i++)
	{
		tmpM = (*this) * tmpM;
	}

	return tmpM;
}


// overload dot-divide and dot-multiply who run just like as they in matlab
template <class Type>
NyxMat<Type> NyxMat<Type>::DotDiv(const NyxMat<Type>& dotDivM)
{
	if ( (m_rows != dotDivM.GetRow()) || (m_cols != dotDivM.GetCol()) )
	{
		std::cerr << "NyxMat::DotDiv: size mismatch." << std::endl;
		exit(3);
	}

	NyxMat<Type> tmpM(m_rows, m_cols);
	for ( int i = 0; i < m_rows; i++)
	{
		for (int j = 0; j < m_cols; j++)
		{
			if (dotDivM(i, j) == 0)
				tmpM(i, j) = 0x7fffffff;
			else
				tmpM(i, j) = m_mat[i][j] / dotDivM(i, j);
		}
	}

	return tmpM;
}



template <class Type>
NyxMat<Type> NyxMat<Type>::DotMulti(const NyxMat<Type>& dotMultiM)
{
	if ( (m_rows != dotMultiM.GetRow()) || (m_cols != dotMultiM.GetCol()) )
	{
		std::cerr << "NyxMat::DotMulti: size mismatch." << std::endl;
		exit(3);
	}

	NyxMat<Type> tmpM(m_rows, m_cols);
	for ( int i = 0; i < m_rows; i++)
	{
		for (int j = 0; j < m_cols; j++)
		{
			tmpM(i, j) = m_mat[i][j] * dotMultiM(i, j);
		}
	}

	return tmpM;
}


// such as A.^2
template <class Type>
NyxMat<Type> NyxMat<Type>::DotPower(int order)
{
	if (order < 0)
	{
		std::cerr << " Only positive integer orders are supported!" << std::endl;
		std::cerr << " Program will quit........." << std::endl;
		exit(3);
	}

	NyxMat <Type> tmpM(m_rows, m_cols);

	for ( int i = 0; i < m_rows; i++)
		for (int j = 0; j < m_cols; j++)
			tmpM(i, j) = static_cast<Type>(pow(static_cast<double>(m_mat[i][j]), order));

	return tmpM;
}



// extend matrix in four directions at the same time (i.e. image dilation)
template <class Type>
NyxMat<Type> NyxMat<Type>::ExtendImg(int OffsetR, int OffsetC)
{
	int i = 0, j = 0;

	if ( (OffsetR < 0) || (OffsetC < 0))
	{
		std::cerr << "NyxMat::ExtendImg: offsets must be positive integers." << std::endl;
		exit(3);
	}
	
	NyxMat<Type> tmpM(m_rows+2*OffsetR, m_cols+2*OffsetC);
	
	// copy original matrix to the center of extended matrix
	for (i = 0; i < m_rows; i++)
	{
		for (j = 0; j < m_cols; j++)
		{
			tmpM(i+OffsetR, j+OffsetC) = m_mat[i][j];
		}
	}
	
	// set values of extended rows
	for (i = 0; i < OffsetR; i++)
	{
		for (j = 0; j < m_cols; j++)
		{
			tmpM(i, j+OffsetC) = m_mat[0][j];
			tmpM(i+m_rows+OffsetR, j+OffsetC) = m_mat[m_rows-1][j];
		}
	}
	
	// set values of extended columns
	for (i = 0; i < m_rows+2*OffsetR; i++)
	{
		for (j = 0; j < OffsetC; j++)
		{
			tmpM(i, j) = tmpM(i, OffsetC);
			tmpM(i, j+OffsetC+m_cols) = tmpM(i, m_cols+OffsetC-1);
		}
	}

	return tmpM;
}


// NyxMat ContractImg(int wOffset, int hOffset);


// minus A --------> -A
template <class Type>
NyxMat<Type> NyxMat<Type>::operator- ()
{
	NyxMat<Type> tmpM(m_rows, m_cols);
	for (int i = 0; i < m_rows; i ++)
		for (int j = 0; j < m_cols; j++)
			tmpM(i, j) = -m_mat[i][j];

	return tmpM;
}


/////////////////////////////////////////////////////////////////////////////////////
// compare A, B -------> a matrix with elements 1 or 0 which represent true or false
template <class Type> template <class Type2>
NyxMat<int> NyxMat<Type>::operator< (const NyxMat<Type2>& compM)
{
	if ( (m_rows != compM.GetRow()) || (m_cols != compM.GetCol()))
	{
		std::cerr << "NyxMat::operator<: size mismatch." << std::endl;
		exit(3);	
	}
	
	NyxMat<int> tmpM(m_rows, m_cols);

	for ( int i = 0; i < m_rows; i++)
	{
		for (int j = 0; j < m_cols; j++)
		{
			if (m_mat[i][j] < compM(i, j))
				tmpM(i, j) = 1;
			else
				tmpM(i, j) = 0;
		}
	}

	return tmpM;
}


template <class Type> template <class Type2>
NyxMat<int> NyxMat<Type>::operator<= (const NyxMat<Type2>& compM)
{
	if ( (m_rows != compM.GetRow()) || (m_cols != compM.GetCol()))
	{
		std::cerr << "NyxMat::operator<: size mismatch." << std::endl;
		exit(3);	
	}
	
	NyxMat<int> tmpM(m_rows, m_cols);

	for ( int i = 0; i < m_rows; i++)
	{
		for (int j = 0; j < m_cols; j++)
		{
			if (m_mat[i][j] <= compM(i, j))
				tmpM(i, j) = 1;
			else
				tmpM(i, j) = 0;
		}
	}

	return tmpM;
}


template <class Type> template <class Type2>
NyxMat<int> NyxMat<Type>::operator< (Type2 cmpNum)
{
	NyxMat<int> tmpM(m_rows, m_cols);

	for ( int i = 0; i < m_rows; i++)
	{
		for (int j = 0; j < m_cols; j++)
		{
			if (m_mat[i][j] < cmpNum)
				tmpM(i, j) = 1;
			else
				tmpM(i, j) = 0;
		}
	}

	return tmpM;
}

template <class Type> template <class Type2>
NyxMat<int> NyxMat<Type>::operator<= (Type2 cmpNum)
{
	NyxMat<int> tmpM(m_rows, m_cols);

	for ( int i = 0; i < m_rows; i++)
	{
		for (int j = 0; j < m_cols; j++)
		{
			if (m_mat[i][j] <= cmpNum)
				tmpM(i, j) = 1;
			else
				tmpM(i, j) = 0;
		}
	}

	return tmpM;
}

///////////////////////////////////////////////////////////////////////////////



template <class Type>
NyxMat<Type>::~NyxMat()
{
	this->DeAllocMem();
}

/*
template <class Type>
template <class Type2> NyxMat<Type>::operator Type2()
{
	Type2 tmpM(rows, cols);
	for (int i = 0; i < rows*cols; i++)
	{
		tmpM(i) = (*this)(i);
	}
	return tmpM;
}
*/

template <class SrcType> template <class DstType>
NyxMat<SrcType>::operator NyxMat<DstType>() const
{
	NyxMat<DstType> tmpM(m_rows, m_cols);
	for (int i = 0; i < m_rows; i++)
		for (int j = 0; j < m_cols; j++)
			tmpM(i, j) = static_cast<DstType>(m_mat[i][j]);

	return tmpM;
}



/****************************************************************************
 *                   Implementions of non-member functions 
 ****************************************************************************/

// overload operator "<<" to output a matrix
template <class Type>
std::ostream& operator<< (std::ostream& out, const NyxMat<Type>& outImg)
{
	for (int i = 0; i < outImg.GetRow(); i++)
	{
		for (int j = 0; j < outImg.GetCol(); j++)
		{
			out << std::setiosflags( std::ios::left)
				<< std::setw(8) << outImg(i, j) << "  ";  // output elements in a readable format
		}
		out << std::endl;
	}
	
	return out;
}

// overload operator ">>" to input a matrix
template <class Type>
std::istream& operator>> (std::istream& in, NyxMat<Type>& inImg)
{
	for (int i = 0; i < inImg.GetRow(); i++)
	{
		std::cout << "Input row " << i << ":" << std::endl;

		for (int j = 0; j < inImg.GetCol(); j++)
		{
			in >> inImg(i, j);                // input elements
		}
	}
	
	return in;
}

// multipying, dividing, adding and substraction between 
// a matrix and a normal number(such as a int...)
// And more, the type of the result matrix is determined 
// by the original matrix, not the number 8-))

template <class Type1, class Type2>
NyxMat<Type1> operator+ (const NyxMat<Type1>& srcM, const NyxMat<Type2>& addM)
{
	int rows = srcM.GetRow();
	int cols = srcM.GetCol();
	if ( (rows != addM.GetRow()) || (cols != addM.GetCol()) )
	{
		std::cerr << "operator+(srcM, addM): size mismatch." << std::endl;
		exit(3);
	}

	NyxMat<Type1> tmpM(rows, cols);  // temporary matrix
	for ( int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			tmpM(i, j) = srcM(i, j) + static_cast<Type1>(addM(i, j));

	return tmpM;
}

template <class Type, class TypeNum>
NyxMat<Type> operator+ (const NyxMat<Type>& srcM, const TypeNum nums)
{
	NyxMat<Type> tmpM(srcM.GetRow(), srcM.GetCol());
	if (tmpM.GetMatrix() == 0)
	{
		tmpM = NyxMat<Type>(srcM.GetRow(), srcM.GetCol());   		// reallocate
		if (tmpM.GetMatrix() == 0)
		{
			std::cerr << "operator+: malloc fail..." << std::endl;
			exit(1);
		}
	}

	for (int i = 0; i < srcM.GetRow(); ++i)
	{
		for (int j = 0; j < srcM.GetCol(); ++j)
		{
			tmpM(i, j) = srcM(i, j) + static_cast<Type>(nums);
		}
	}
	return tmpM;
}


template <class Type1, class Type2>
NyxMat<Type1> operator- (const NyxMat<Type1>& srcM, const NyxMat<Type2>& subM)
{
	int rows = srcM.GetRow();
	int cols = srcM.GetCol();
	if ( (rows != subM.GetRow()) || (cols != subM.GetCol()) )
	{
		std::cerr << "operator-(srcM, subM): size mismatch." << std::endl;
		exit(3);
	}

	NyxMat<Type1> tmpM(rows, cols);  // temporary matrix
	for ( int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			tmpM(i, j) = srcM(i, j) - static_cast<Type1>(subM(i, j));

	return tmpM;
}


template <class Type, class TypeNum> 
NyxMat<Type> operator- (const NyxMat<Type>& srcM, const TypeNum nums)
{
	NyxMat<Type> tmpM(srcM.GetRow(), srcM.GetCol());
	if (tmpM.GetMatrix() == 0)
	{
		tmpM = NyxMat<Type>(srcM.GetRow(), srcM.GetCol());   		// reallocate
		if (tmpM.GetMatrix() == 0)
		{
			std::cerr << "operator-: tmpM malloc fail." << std::endl;
			exit(1);
		}
	}

	for (int i = 0; i < srcM.GetRow(); ++i)
	{
		for (int j = 0; j < srcM.GetCol(); ++j)
		{
			tmpM(i, j) = srcM(i, j) - static_cast<Type>(nums);
		}
	}
	return tmpM;
}


template <class Type1, class Type2>
NyxMat<Type1> operator* (const NyxMat<Type1>& srcM, const NyxMat<Type2>& multiM)
{
	int rows = srcM.GetRow();
	int cols = srcM.GetCol();

	if( cols != multiM.GetRow() ) 
	{
		std::cerr << "operator*(srcM, multiM): size mismatch." << std::endl;
		exit(3);
	}
	
    NyxMat<Type1> tmpM(rows, multiM.GetCol());

    for(int i = 0; i < rows; i++)
	{
		for(int j = 0; j < multiM.GetCol(); j++)
        {
            tmpM(i, j) = 0;
            for(int k = 0; k < cols; k++)
				tmpM(i, j) += srcM(i, k) * static_cast<Type1>(multiM(k, j));
        }
	}

    return tmpM;
}


template <class Type, class TypeNum> 
NyxMat<Type> operator* (const NyxMat<Type>& srcM, const TypeNum nums)
{
	NyxMat<Type> tmpM(srcM.GetRow(), srcM.GetCol());
	if (tmpM.GetMatrix() == 0)
	{
		tmpM = NyxMat<Type>(srcM.GetRow(), srcM.GetCol());   		// reallocate
		if (tmpM.GetMatrix() == 0)
		{
			std::cerr << "operator-: tmpM malloc fail." << std::endl;
			exit(1);
		}
	}

	for (int i = 0; i < srcM.GetRow(); ++i)
	{
		for (int j = 0; j < srcM.GetCol(); ++j)
		{
			tmpM(i, j) = srcM(i, j) * static_cast<Type>(nums);
		}
	}

	return tmpM;
}

//template <class Type, class TypeNum> 
//NyxMat<Type> operator* (TypeNum nums, NyxMat<Type>& srcM)
//{
//	return srcM*nums;
//}


// matrix divition is very hard... leave space for supplement here
//template <class Type1, class Type2>
//NyxMat<Type1> operator/ (const NyxMat<Type1>& srcM, const NyxMat<Type2>& divM)
//{
//}


template <class Type, class TypeNum> 
NyxMat<Type> operator/ (const NyxMat<Type>& srcM, const TypeNum nums)
{
	if (nums == 0)
		throw std::out_of_range("Divided by 0");
	
	NyxMat<Type> tmpM(srcM.GetRow(), srcM.GetCol());
	if (tmpM.GetMatrix() == 0)
	{
		tmpM = NyxMat<Type>(srcM.GetRow(), srcM.GetCol());   		// reallocate
		if (tmpM.GetMatrix() == 0)
		{
			std::cerr << "operator/: tmpM Malloc fail." << std::endl;
			exit(1);
		}
	}

	for (int i = 0; i < srcM.GetRow(); ++i)
	{
		for (int j = 0; j < srcM.GetCol(); ++j)
		{
			tmpM(i, j) = srcM(i, j) / static_cast<Type>(nums);
		}
	}
	return tmpM;
}


// exp() always returns NyxMat<double>
template <class Type>
NyxMat<double> MatExp(const NyxMat<Type>& srcM)
{
	NyxMat<double> tmpM(srcM.GetRow(), srcM.GetCol());
	if (tmpM.GetMatrix() == 0)
	{
		tmpM = NyxMat<Type>(srcM.GetRow(), srcM.GetCol());   		// reallocate
		if (tmpM.GetMatrix() == 0)
		{
			std::cerr << "exp: tmpM Malloc fail." << std::endl;
			exit(1);
		}
	}

	for (int i = 0; i < srcM.GetRow(); ++i)
	{
		for (int j = 0; j < srcM.GetCol(); ++j)
		{
			tmpM(i, j) = exp(static_cast<double>(srcM(i, j)));
		}
	}
	return tmpM;
}


// cross multiplication: only for 2D or 3D vectors
template <class Type>
NyxMat<Type> CrossMulti(const NyxMat<Type>& vec1, const NyxMat<Type>& vec2)
{
	int row1 = vec1.GetRow();
	int col1 = vec1.GetCol();
	int row2 = vec2.GetRow();
	int col2 = vec2.GetCol();

	Type a;

	if (row1==1 && row2==1 && col1==col2)
	{
		if (col1==2)
		{
			a = vec1(0,0)*vec2(0,1)-vec2(0,0)*vec1(0,1);
			if ( a>=0 && a<=powf(10, -9))
				a = powf(10, -9);
			else if (a<0 && a>-powf(10, -9))
				a = -powf(10, -9);

			NyxMat<Type> tmpM(1, 2);
			
			tmpM(0, 0) = (vec1(0,1)-vec2(0,1))/a;
			tmpM(0, 1) = (vec2(0,0)-vec1(0,0))/a;
			return tmpM;
		}
		else if (col1==3)
		{
			NyxMat<Type> tmpM(1, 3);
			tmpM(0, 0) = vec1(0,1)*vec2(0,2)-vec1(0,2)*vec2(0,1);
			tmpM(0, 1) = vec1(0,2)*vec2(0,0)-vec2(0,2)*vec1(0,0);
			tmpM(0, 2) = vec1(0,0)*vec2(0,1)-vec2(0,0)*vec1(0,1);
			return tmpM;
		}
		else
		{
		  std::cerr << "crossMulti error: " << std::endl <<
		    "	only vectors of 1x2, 2x1, 1x3, 3x1 are supported temporarily." << std::endl;
			exit(-1);
		}
	}
	else if (col1==1 && col2==1 && row1==row2)
	{
		if (row1==2)
		{
			a = vec1(0,0)*vec2(1,0)-vec2(0,0)*vec1(1,0);
			if ( a>=0 && a<=powf(10, -9))
				a = powf(10, -9);
			else if (a<0 && a>-powf(10, -9))
				a = -powf(10, -9);

			NyxMat<Type> tmpM(2, 1);
			
			tmpM(0, 0) = (vec1(1,0)-vec2(1,0))/a;
			tmpM(1, 0) = (vec2(0,0)-vec1(0,0))/a;
			return tmpM;
		}
		else if (row1==3)
		{
			NyxMat<Type> tmpM(3, 1);
			tmpM(0, 0) = vec1(1,0)*vec2(2,0)-vec1(2,0)*vec2(1,0);
			tmpM(1, 0) = vec1(2,0)*vec2(0,0)-vec2(2,0)*vec1(0,0);
			tmpM(2, 0) = vec1(0,0)*vec2(1,0)-vec2(0,0)*vec1(1,0);
			return tmpM;
		}
		else
		{
		  std::cerr << "crossMulti error: " << std::endl <<
		    "	only vectors of 1x2, 2x1, 1x3, 3x1 are supported temporarily." << std::endl;
			exit(-1);
		}
	}
	else
	{
	  std::cerr << "crossMulti error: " << std::endl <<
	    "	Dimensions of the two vectors must be same." << std::endl;
		exit(-1);
	}
	return vec1;
}


// sum of all the elements in the matrix
template <class Type>
Type MatSum(const NyxMat<Type>& srcM)
{
	Type rlt = 0;
	
	for (int i = 0; i < srcM.GetRow(); ++i)
	{
		for (int j = 0; j < srcM.GetCol(); ++j)
		{
			rlt += srcM(i, j);
		}
	}

	return rlt;
}


// 2-norm of matrix
template <class Type>
double Mat2Norm(const NyxMat<Type>& srcM)
{
	double sum = 0;   // "Type" is converted into "double" implicitly

	for (int i = 0; i < srcM.GetRow(); ++i)
	{
		for (int j = 0; j < srcM.GetCol(); ++j)
		{
			sum += srcM(i, j) * srcM(i, j);
		}
	}
	return sqrt(sum);
}


// two-dimension mediate Filters with window size of wndSize x wndSize
template <class Type>
NyxMat<Type> MedFilt2(NyxMat<Type>& srcM, int wndSize)
{
	int offset = wndSize / 2;
	NyxMat<Type> nImg(srcM.ExtendImg(offset, offset));
	NyxMat<Type> rlt(srcM.GetRow(), srcM.GetCol(), 5);

	// a small window used in mediate filters
	NyxMat<Type> patch(wndSize, wndSize);
	
	for ( int i = 0; i < srcM.GetRow(); i++)
	{
		for ( int j = 0; j < srcM.GetCol(); j++)
		{
			// sort the patch window descendly using insert-sortion
			patch = NyxMat<Type>(nImg, i, wndSize, j, wndSize);  // crop a patch
			InsertSort(patch);

			rlt(i, j) = patch(wndSize*wndSize/2);	
		}
	}

	return rlt;
}


template <class Type>
void InsertSort(NyxMat<Type>& patch)
{
	int wndSize = patch.GetCol();

	for (int k = 0; k < wndSize*wndSize-1; k++)
	{
		int point = k;
		for (int m = k+1; m < wndSize*wndSize; m++)
		{
			if (patch(m) > patch(point))
			{
				point = m;
			}
		}
		if ( point != k)
		{
			Type tmp = patch(k);
			patch(k) = patch(point);
			patch(point) = tmp;
		}
	}
}
/*
template <class Type>
Zeros<Type>::Zeros(unsigned int rows, unsigned int cols)
	: NyxMat<Type>(rows, cols, 0)
{
}

template <class Type>
Ones<Type>::Ones(unsigned int rows, unsigned int cols)
	: NyxMat<Type>(rows, cols, 1)
{
}

template <class Type>
Diag<Type>::Diag(unsigned int rows, unsigned int cols, Type nums)
	: NyxMat<Type>(rows, cols, nums, true)
{
}
*/


// matrix inversioon
template <class Type>
NyxMat<double> MatInv(const NyxMat<Type> &A)
{
	if (A.GetRow() != A.GetCol())
	{
	  std::cerr << "MatInv error: Only square matrix are supported." << std::endl;
		exit(-1);
	}
	int order = A.GetRow();

    // get the determinant of matrix A
    double det = 1.0/MatDet(A);

	NyxMat<double> Y(order, order);   // store the result

    // memory allocation
	NyxMat<Type> nxminor(order-1, order-1);

    for(int j = 0; j < order; ++j)
    {
        for(int i = 0; i < order; ++i)
        {
            // get the co-factor (matrix) of A(j,i)
            GetMinor(A, nxminor, j, i);
            Y(i, j) = MatDet(nxminor) * det;

            if( (i + j)%2 == 1 )
                Y(i, j) = -Y(i, j);
        }
    }
	
	return Y;
}



// calculate the cofactor of element (row,col)
template <class Type>
void GetMinor(const NyxMat<Type> &src, NyxMat<Type> &dest, int row, int col)
{
	int order = src.GetRow();

    // indicate which col and row is being copied to dest
    int colCount = 0, rowCount = 0;

    for(int i = 0; i < order; ++i)
    {
        if( i != row )
        {
            colCount = 0;
            for(int j = 0; j < order; ++j)
            {
                // when j is not the element
                if( j != col )
                {
                    dest(rowCount, colCount) = src(i, j);
                    colCount++;
                }
            }
            rowCount++;
        }
    }
}


// Calculate the determinant recursively.
template <class Type>
Type MatDet(const NyxMat<Type> &mat)
{
	if (mat.GetRow() != mat.GetCol())
	{
		std::cerr << "MatDet error: matrix must be square." << std::endl;
		exit(-1);
	}

	int order = mat.GetRow();

    // order must be >= 0
	// stop the recursion when matrix is a single element
    if( order == 1 )
        return mat(0, 0);

    // the determinant value
    Type det = 0;

    // allocate the cofactor matrix
	NyxMat<Type> nxminor(order-1, order-1);

    for(int i = 0; i < order; i++ )
    {
        // get nxminor of element (0,i)
        GetMinor(mat, nxminor, 0, i);

        // the recusion is here!
        det += pow(-1.0, i) * mat(0, i) * MatDet(nxminor);
    }

    return det;
}




#endif
