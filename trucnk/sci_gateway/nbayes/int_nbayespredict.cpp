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

#include "naivebayes.h"
#include "mattransform.h"


extern "C"
{
	#include "stack-c.h"
	#include "Scierror.h"
	#include "sciprint.h"
	
	#include <stdio.h>
	
	#undef CreateVar
	#define CreateVar(n,ct,mx,nx,lx) if(! C2F(createvar)((c_local=n,&c_local),ct,mx,nx,lx, 1L))\
	        { return 0;  }
	        
	int int_nbayespredict(char* fname)
	{	
		int m, n, l;
		const size_t mat_size = sizeof(CvMat*);
		int nclasses;
		
		char* pStr;
		SciIntMat Mat;
		
		CvMat* samples;	 //data type of CV_32FC1
		CvMat* responses;  //data type of CV_32FC1
	
	    int  var_count, var_all;
	    CvMat*  var_idx;  //data type of CV_32SC1	
	    CvMat*  cls_labels;  //data type of CV_32SC1
	    CvMat** count;  //data type of CV_32SC1
	    CvMat** sum;
	    CvMat** productsum;
	    CvMat** avg;
	    CvMat** inv_eigen_values;
	    CvMat** cov_rotate_mats;
	    CvMat*  c;
	    
	    CvNormalBayesClassifier nbc;
	    
		if(Rhs == 0)
		{
			sciprint("Usage:\ntest_labels = nbayespredict(test_data, model)\n");
			return -1;
		}
		
		CheckRhs(2, 2);
		CheckLhs(1, 1);
		
		//prepare input data for prediction
		//test_data
		samples = SciMat2CvMat(1, 5);
		responses = cvCreateMat(1, samples->rows, CV_32FC1);
		cvSetZero(responses);
		
		//read mlist
		GetListRhsVar(2, 1, MATRIX_OF_STRING_DATATYPE, &m, &n, &l);	
		pStr = cstk(l);
		//var_count
		GetListRhsVar(2, 2, MATRIX_OF_VARIABLE_SIZE_INTEGER_DATATYPE, &m, &n, &Mat);
		var_count = ((int*)Mat.D)[0];
		//var_all 
		GetListRhsVar(2, 3, MATRIX_OF_VARIABLE_SIZE_INTEGER_DATATYPE, &m, &n, &Mat);
		var_all = ((int*)Mat.D)[0];
		//var_idx
		GetListRhsVar(2, 4, MATRIX_OF_VARIABLE_SIZE_INTEGER_DATATYPE, &m, &n, &Mat);
		if((m != 0) && (n != 0))  //if var_idx is not null
		{
			var_idx = cvCreateMat(m, n, CV_32SC1);
			CvMat* var_idx_tmp = cvCreateMat(n, m, CV_32SC1);
			var_idx_tmp->data.i = (int*)Mat.D;
			cvTranspose(var_idx_tmp, var_idx);
			cvReleaseMat(&var_idx_tmp);		
		}
		else  //var_idx is null
		{
			var_idx = NULL;
		}
		//cls_labels
		GetListRhsVar(2, 5, MATRIX_OF_VARIABLE_SIZE_INTEGER_DATATYPE, &m, &n, &Mat);
		cls_labels = cvCreateMat(m, n, CV_32SC1);
		nclasses = cls_labels->cols;
		CvMat* cls_labels_tmp = cvCreateMat(n, m, CV_32SC1);
		cls_labels_tmp->data.i = (int*)Mat.D;
		cvTranspose(cls_labels_tmp, cls_labels);
		cvReleaseMat(&cls_labels_tmp);
		//count
		GetListRhsVar(2, 6, MATRIX_OF_VARIABLE_SIZE_INTEGER_DATATYPE, &m, &n, &Mat);
		CvMat* tmp1 = cvCreateMat(m, n, CV_32SC1);
		tmp1->data.i = (int*)Mat.D;		
		count = Array2CvMat(tmp1, nclasses);
		cvReleaseMat(&tmp1);
		//sum
		GetListRhsVar(2, 7, MATRIX_OF_DOUBLE_DATATYPE, &m, &n, &l);
		CvMat* tmp2 = cvCreateMat(m, n, CV_64FC1);
		tmp2->data.db = stk(l);
		sum = Array2CvMat(tmp2, nclasses);
		cvReleaseMat(&tmp2);
		//productsum
		GetListRhsVar(2, 8, MATRIX_OF_DOUBLE_DATATYPE, &m, &n, &l);
		CvMat* tmp3 = cvCreateMat(m, n, CV_64FC1);
		tmp3->data.db = stk(l);
		productsum = Array2CvMat(tmp3, nclasses);
		cvReleaseMat(&tmp3);
		//avg
		GetListRhsVar(2, 9, MATRIX_OF_DOUBLE_DATATYPE, &m, &n, &l);
		//printf("m-avg: %d n-avg: %d\n", m, n);						
		CvMat* tmp4 = cvCreateMat(m, n, CV_64FC1);
		tmp4->data.db = stk(l);
		avg = Array2CvMat(tmp4, nclasses);
		cvReleaseMat(&tmp4);	
		//inv_eigen_values
		GetListRhsVar(2, 10, MATRIX_OF_DOUBLE_DATATYPE, &m, &n, &l);
		CvMat* tmp5 = cvCreateMat(m, n, CV_64FC1);
		tmp5->data.db = stk(l);
		inv_eigen_values = Array2CvMat(tmp5, nclasses);
		cvReleaseMat(&tmp5);		
		//cov_rotate_mats
		GetListRhsVar(2, 11, MATRIX_OF_DOUBLE_DATATYPE, &m, &n, &l);
		CvMat* tmp6 = cvCreateMat(m, n, CV_64FC1);
		tmp6->data.db = stk(l);
		cov_rotate_mats = Array2CvMat(tmp6, nclasses);
		cvReleaseMat(&tmp6);	
		//c
		GetListRhsVar(2, 12, MATRIX_OF_DOUBLE_DATATYPE, &m, &n, &l);
		//printf("m-c: %d n-c: %d\n", m, n);								
		CvMat* c_tmp = cvCreateMat(n, m, CV_64FC1);
		c_tmp->data.db = stk(l);
		c = cvCreateMat(m, n, CV_64FC1);
		cvTranspose(c_tmp, c);
		cvReleaseMat(&c_tmp);	
		
		//normal bayes model
		nbc.var_count = var_count;
		nbc.var_all = var_all;
		nbc.var_idx = var_idx;
		nbc.cls_labels = cls_labels;
		nbc.count = count;
		nbc.sum = sum;
		nbc.productsum = productsum;
		nbc.avg = avg;
		nbc.inv_eigen_values = inv_eigen_values;
		nbc.cov_rotate_mats = cov_rotate_mats;
		nbc.c = c;
		
		//predict 
		nbc.predict(samples, responses);
		
		CvMat2SciMat(responses, Rhs+1);
		
		LhsVar(1) = Rhs+1;
		
		//release memory
		cvReleaseMat(&samples);
		cvReleaseMat(&responses);
		
		if(cls_labels)
		{
	        for(int cls = 0; cls < cls_labels->cols; cls++)
	        {
	            cvReleaseMat(&count[cls] );
	            cvReleaseMat(&sum[cls] );
	            cvReleaseMat(&productsum[cls] );
	            cvReleaseMat(&avg[cls] );
	            cvReleaseMat(&inv_eigen_values[cls] );
	            cvReleaseMat(&cov_rotate_mats[cls] );
	        }
        }
        
/*    	if(cls_labels)
	    	cvReleaseMat( &cls_labels );
	    if(var_idx)
	    	cvReleaseMat( &var_idx );
	    cvFree(&count);
	    cvFree(&sum);
	    cvFree(&productsum);
	    cvFree(&avg);
	    cvFree(inv_eigen_values);
	    cvFree(cov_rotate_mats);	
	    cvReleaseMat( &c );*/					
			
		return 0;
	}
	
}	

