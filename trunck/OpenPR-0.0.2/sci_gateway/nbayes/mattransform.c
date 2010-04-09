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

#include "mattransform.h"
#include "common.h"


/***********************************************************
* convert SCI matrix to OpenCV CvMat 
* change the data order from column-wise to row-wise
***********************************************************/ 
CvMat* SciMat2CvMat(int nPos, int flag)
{
	CvMat* pMat = NULL;
	int m, n, l;
    int i, j, k;
    
    GetRhsVar(nPos, MATRIX_OF_DOUBLE_DATATYPE, &m, &n, &l);
    //double* pSrc = stk(l);
    
    //in cxtypes.h
    // CV_8U   0
	// CV_8S   1
	// CV_16U  2
	// CV_16S  3
	// CV_32S  4
	// CV_32F  5
	// CV_64F  6
	// CV_USRTYPE1 7
    
    switch(flag)
    {
    case 0:
    	pMat = cvCreateMat(m, m, CV_8UC1);
	    k = 0;
	    for(i = 0; i < n; i++)
	    {
	    	for(j = 0; j < m; j++)
	    	{
	    		pMat->data.ptr[j*n+i] = (unsigned char)stk(l)[k++];//(unsigned char)pSrc[k++];
	    	}
	    }
	    break;
    case 4:
    	pMat = cvCreateMat(m, m, CV_32SC1);
	    k = 0;
	    for(i = 0; i < n; i++)
	    {
	    	for(j = 0; j < m; j++)
	    	{
	    		pMat->data.i[j*n+i] = (int)stk(l)[k++];//(int)pSrc[k++];
	    	}
	    }
	    break;
    case 5:
	    pMat = cvCreateMat(m, n, CV_32FC1);  
	    k = 0;
	    for(i = 0; i < n; i++)
	    {
	    	for(j = 0; j < m; j++)
	    	{
	    		pMat->data.fl[j*n+i] = (float)stk(l)[k++];//(float)pSrc[k++];
	    	}
	    }
	    break;
	default:
		break;
    }
        
    return pMat;
}

/***********************************************************
* convert OpenCV CvMat (CV_32FC1 or CV_32SC1) to SCI matrix
* change the data order from row-wise to column-wise
***********************************************************/
void CvMat2SciMat(CvMat* pMat, int nPos)
{
	int rows = pMat->rows;
	int cols = pMat->cols;
	
	if(CV_MAT_TYPE(pMat->type) == CV_32FC1)
	{
		float *pData = pMat->data.fl;
		Create2DFloatMat(nPos, rows, cols, pData);		
	}
	
	if(CV_MAT_TYPE(pMat->type) == CV_32SC1)
	{
		int *pData = pMat->data.i;
		int nType = I_INT32;
		Create2DIntMat(nPos, rows, cols, pData, nType);
	}	
}


void CvMat2Array(CvMat** pMat, void* pArray, int dim)
{
	int r, cl, dm;
	
	int i, j, d;
	int k;
	
	unsigned char* p_uchar = NULL;
	int *p_int = NULL;
	double* p_db = NULL;
	
	r = pMat[0]->rows;
	cl = pMat[0]->cols;
	dm = dim;	
	
	switch(CV_MAT_TYPE(pMat[0]->type))
	{
	case CV_8UC1:
		p_uchar = (unsigned char*)pArray;
		k = 0;	
		for(d = 0; d < dm; d++)
		{
			for(j = 0; j < cl; j++)
			{
				for(i = 0; i < r; i++)
				{
					p_uchar[k++] = CV_MAT_ELEM((*(pMat[d])), unsigned char, i, j);
				}
			}
		}
		break;
	case CV_32SC1:
		p_int = (int*)pArray;
		k = 0;	
		for(d = 0; d < dm; d++)
		{
			for(j = 0; j < cl; j++)
			{
				for(i = 0; i < r; i++)
				{
					p_int[k++] = CV_MAT_ELEM((*(pMat[d])), int, i, j);
				}
			}
		}
		break;
	case CV_64FC1:
		p_db = (double*)pArray;
		k = 0;	
		for(d = 0; d < dm; d++)
		{
			for(j = 0; j < cl; j++)
			{
				for(i = 0; i < r; i++)
				{
					p_db[k++] = CV_MAT_ELEM((*(pMat[d])), double, i, j);
				}
			}
		}
		break;
	default:
		break;
	}
}

CvMat** Array2CvMat(CvMat* pArray, int nclasses)
{
	//type:
	//CV_32SC1 4
	//CV_64FC1 6
	
	CvMat ** ppMat = (CvMat**)cvAlloc(nclasses*sizeof(CvMat*));
	
	int* pData_int;
	double* pData_db;
	
	int i, j, k;
	int idx;
	
	int m, n, dim;
	int rows;
	m = pArray->rows;
	n = pArray->cols;
	dim = (int)(m/n);

	if(dim == nclasses)
	{
		rows = n;
	}
	else
	{
		rows = 1;
	}		

	for(k = 0; k < nclasses; k++)
	{
		ppMat[k] = cvCreateMat(rows, n, CV_MAT_TYPE(pArray->type));
	}		
	
	switch(CV_MAT_TYPE(pArray->type))
	{
	case CV_32SC1:
		pData_int = pArray->data.i;	
		idx = 0;	
		for(j = 0; j < n; j++)
		{
			for(k = 0; k < nclasses; k++)
			{
				for(i = 0; i < rows; i++)
				{
					CV_MAT_ELEM((*ppMat[k]), int, i, j) = pData_int[idx++];
				}
			}
		}
		break;
	case CV_64FC1:
		pData_db = pArray->data.db;
		idx = 0;	
		for(j = 0; j < n; j++)
		{
			for(k = 0; k < nclasses; k++)
			{
				for(i = 0; i < rows; i++)
				{
					CV_MAT_ELEM((*ppMat[k]), double, i, j) = pData_db[idx++];
				}
			}
		}
		break;
	default:
		break;
	}
	
	return ppMat;	
}