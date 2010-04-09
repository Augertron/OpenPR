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
	#include "sciprint.h"
	#include "Scierror.h"
	
	#include <stdio.h>
	
	#undef CreateVar
	#define CreateVar(n,ct,mx,nx,lx) if(! C2F(createvar)((c_local=n,&c_local),ct,mx,nx,lx, 1L))\
	        { return 0;  }
	        
	int int_nbayestrain(char* fname)
	{
		int m, n, l, l2;
		int i, j, d;
		int k;
		
		CvMat* train_data = NULL;  //data type of CV_32FC1
		CvMat* train_labels = NULL;  //data type of CV_32SC1 or CV_32FC1 (here use CV_32FC1)
		CvMat* var_idx = NULL;  //data type of CV_32SC1 or CV_8UC1 (here use CV_8UC1)
		CvMat* sample_idx = NULL;  //data type of CV_32SC1 or CV_8UC1 (here use CV_8UC1)
		int update = 0;

		CvNormalBayesClassifier nbc;
		
		SciIntMat mat;
		
		//arrays to store the model
//		double* var_count = NULL;
//		double* var_all = NULL;
		int* var_count = NULL;
		int* var_all = NULL;
		int* var_idx_int = NULL;    //the type is CV_32SC1
		int* cls_labels_int = NULL;   //the type is CV_32SC1
		int* count = NULL;
		double* sum = NULL;
		double* productsum = NULL;
		double* avg = NULL;
		double* inv_eigen_values = NULL;
		double* cov_rotate_mats = NULL;
		double* c = NULL;
		
		char* Str[]={"nbc_model","var_count","var_all","var_idx","cls_labels","count","sum","productsum","avg","inv_eigen_values","cov_rotate_mats","c"};
		
		if(Rhs == 0)
		{
			sciprint("Usage:\nmodel = nbayestrain(train_data, train_labels[, var_idx[, sample_idx]])\n");
			return -1;
		}
		
		CheckRhs(2, 4);
		CheckLhs(1, 1);
		
		//prepare the input data
		train_data = SciMat2CvMat(1, 5);
		train_labels = SciMat2CvMat(2, 5);	  //use data type CV_32FC1

		switch(Rhs)
		{
		case 3:
			var_idx = SciMat2CvMat(3, 0);  //use data type CV_8UC1
			if((var_idx->rows != 1) && (var_idx->cols !=1))
			{
				Scierror(999, "%s: var_idx mask should be 1-dimensional\n", fname);
				return -1;
			}
			if((var_idx->rows+var_idx->cols-1) != train_data->cols)
			{
				Scierror(999, "%s: var_idx mask should contain as many elements as the total number of input data features\n", fname);
				return -1;
			}
			break;
		case 4:
			var_idx = SciMat2CvMat(3, 0);  //use data type CV_8UC1
			sample_idx = SciMat2CvMat(4, 0);  //use data type CV_8UC1
			if((var_idx->rows != 1) && (var_idx->cols !=1))
			{
				Scierror(999, "%s: var_idx mask should be 1-dimensional\n", fname);
				return -1;
			}
			if((var_idx->rows+var_idx->cols-1) != train_data->cols)
			{
				Scierror(999, "%s: var_idx mask should contain as many elements as the total number of input data features\n", fname);
				return -1;
			}

			if((sample_idx->rows != 1) && (sample_idx->cols !=1))
			{
				Scierror(999, "%s: sample_idx mask should be 1-dimensional\n", fname);
				return -1;
			}
			if((sample_idx->rows+sample_idx->cols-1) != train_data->rows)
			{
				Scierror(999, "%s: sample_idx mask should contain as many elements as the total number of input data samples\n", fname);
				return -1;
			}					
			break;
/*		case 5:
			var_idx = SciMat2CvMat(3, 0); 
			sample_idx = SciMat2CvMat(4, 0);
			GetRhsVar(5, MATRIX_OF_DOUBLE_DATATYPE, &m, &n, &l);
			update = (int)(*stk(l));
			break;*/
		default:
			break;
		}
		
		//training
		nbc.train(train_data, train_labels, var_idx, sample_idx, update);
		
		//the model is to be stored in mlist data structure
//		var_count = new double;
//		var_all = new double;
		var_count = new int;
		var_all = new int;
		var_count[0] = nbc.var_count;
		var_all[0] = nbc.var_all;

		count = new int[(nbc.cls_labels->cols)*(nbc.var_count)];								//should be stored as hypermat
		sum = new double[(nbc.cls_labels->cols)*(nbc.var_count)];								//should be stored as hypermat
		productsum = new double[(nbc.cls_labels->cols)*(nbc.var_count)*(nbc.var_count)];    	//should be stored as hypermat
		avg = new double[(nbc.cls_labels->cols)*(nbc.var_count)];								//should be stored as hypermat
		inv_eigen_values = new double[(nbc.cls_labels->cols)*(nbc.var_count)];					//should be stored as hypermat
		cov_rotate_mats = new double[(nbc.cls_labels->cols)*(nbc.var_count)*(nbc.var_count)];  	//should be stored as hypermat
		c = new double[nbc.cls_labels->cols];
		
		//create mlist and store data of the model	
		m = 12;
		n = 1;
		CreateVar(Rhs+1, MATRIX_ORIENTED_TYPED_LIST_DATATYPE, &m, &n, &l);
		m = 1;
		n = 12;
		CreateListVarFromPtr(Rhs+1, 1, MATRIX_OF_STRING_DATATYPE, &m, &n, Str);
		
		mat.m = 1;
		mat.n = 1;
		mat.l = -1;
		mat.it = I_INT32;
		mat.D = var_count;
		CreateListVarFromPtr(Rhs+1, 2, MATRIX_OF_VARIABLE_SIZE_INTEGER_DATATYPE, &(mat.m), &(mat.n), &mat);
		mat.D = var_all;
		CreateListVarFromPtr(Rhs+1, 3, MATRIX_OF_VARIABLE_SIZE_INTEGER_DATATYPE, &(mat.m), &(mat.n), &mat);		

		if(nbc.var_idx)
		{
			var_idx_int = new int[(nbc.var_idx->rows)*(nbc.var_idx->cols)];
			CvMat2Array(&(nbc.var_idx), var_idx_int, 1);

			mat.m = nbc.var_idx->rows;
			mat.n = nbc.var_idx->cols;
			mat.l = -1;
			mat.it = I_INT32;
			mat.D = var_idx_int;

			CreateListVarFromPtr(Rhs+1, 4, MATRIX_OF_VARIABLE_SIZE_INTEGER_DATATYPE, &(mat.m), &(mat.n), &mat);
		}
		else
		{
			mat.m = 0;
			mat.n = 0;
			mat.l = -1;
			mat.it = I_INT32;
			mat.D = 0;
			CreateListVarFromPtr(Rhs+1, 4, MATRIX_OF_VARIABLE_SIZE_INTEGER_DATATYPE, &(mat.m), &(mat.n), &mat);
		}
		
		if(nbc.cls_labels)
		{
			cls_labels_int = new int[(nbc.cls_labels->rows)*(nbc.cls_labels->cols)];
			CvMat2Array(&(nbc.cls_labels), cls_labels_int, 1);
			
			mat.m = nbc.cls_labels->rows;
			mat.n = nbc.cls_labels->cols;
			mat.l = -1;
			mat.it = I_INT32;
			mat.D = cls_labels_int;

			CreateListVarFromPtr(Rhs+1, 5, MATRIX_OF_VARIABLE_SIZE_INTEGER_DATATYPE, &(mat.m), &(mat.n), &mat);
		}
		else
		{
			mat.m = 0;
			mat.n = 0;
			mat.l = -1;
			mat.it = I_INT32;
			mat.D = 0;
			CreateListVarFromPtr(Rhs+1, 5, MATRIX_OF_VARIABLE_SIZE_INTEGER_DATATYPE, &(mat.m), &(mat.n), &mat);
		}
	
		CvMat2Array(nbc.count, count, nbc.cls_labels->cols);
		mat.m = nbc.cls_labels->cols;
		mat.n = nbc.var_count;
		mat.l = -1;
		mat.it = I_INT32;
		mat.D = count;
		CreateListVarFromPtr(Rhs+1, 6, MATRIX_OF_VARIABLE_SIZE_INTEGER_DATATYPE, &(mat.m), &(mat.n), &mat);
		
		CvMat2Array(nbc.sum, sum, nbc.cls_labels->cols);
		m = nbc.cls_labels->cols;
		n = nbc.var_count;
		CreateListVarFromPtr(Rhs+1, 7, MATRIX_OF_DOUBLE_DATATYPE, &m, &n, &sum);
		
		CvMat2Array(nbc.productsum, productsum, nbc.cls_labels->cols);
		m = (nbc.cls_labels->cols)*(nbc.var_count);
		n = nbc.var_count;
		CreateListVarFromPtr(Rhs+1, 8, MATRIX_OF_DOUBLE_DATATYPE, &m, &n, &productsum);
		
		CvMat2Array(nbc.avg, avg, nbc.cls_labels->cols);
		m = nbc.cls_labels->cols;
		n = nbc.var_count;
		CreateListVarFromPtr(Rhs+1, 9, MATRIX_OF_DOUBLE_DATATYPE, &m, &n, &avg);
		
		CvMat2Array(nbc.inv_eigen_values, inv_eigen_values, nbc.cls_labels->cols);
		m = nbc.cls_labels->cols;
		n = nbc.var_count;
		CreateListVarFromPtr(Rhs+1, 10, MATRIX_OF_DOUBLE_DATATYPE, &m, &n, &inv_eigen_values);
	
		CvMat2Array(nbc.cov_rotate_mats, cov_rotate_mats, nbc.cls_labels->cols);	
		m = (nbc.cls_labels->cols)*(nbc.var_count);
		n = nbc.var_count;
		CreateListVarFromPtr(Rhs+1, 11, MATRIX_OF_DOUBLE_DATATYPE, &m, &n, &cov_rotate_mats);
		
		CvMat2Array(&(nbc.c), c, 1);	
		m = 1;
		n = nbc.cls_labels->cols;
		CreateListVarFromPtr(Rhs+1, 12, MATRIX_OF_DOUBLE_DATATYPE, &m, &n, &c);
		
		LhsVar(1) = Rhs+1;
		
		//release memory
		cvReleaseMat(&train_data);
		cvReleaseMat(&train_labels);
		
		if(var_idx)
		{
			cvReleaseMat(&var_idx);
		}
		if(sample_idx)
		{
			cvReleaseMat(&sample_idx);
		}	
		if(var_idx_int)
		{
			delete []var_idx_int;
		}
		if(cls_labels_int)
		{
			delete []cls_labels_int;
		}
		
		delete []count;
		delete []sum;
		delete []productsum;
		delete []avg;
		delete []inv_eigen_values;
		delete []cov_rotate_mats;
		delete []c;
		
		return 0;		
	}
	
}
