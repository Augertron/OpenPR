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

#include "ml.h"

#include "common.h"


extern "C"{


	#include "sciprint.h"
	#include "Scierror.h"


	int int_knearest(char *fname)
	{
		if(Rhs == 0)
		{
			sciprint("Usage:\ntest_labels = knearest(train_labels, train_data, test_data, k[, is_regression])\n");
			return -1;
		}

		CheckRhs(4, 5);
		CheckLhs(1, 1);
		
		//_train_data and _test_data should be of CV_32FC1 type
		//_train_labels and _test_labels should be of CV_32FC1 type of CV_32SC1(in case of classification)	
		if((VarType(1) != 1 && VarType(1) != 8) || VarType(2) != 1 || VarType(3) != 1)
		{
			Scierror(999, "%s: argument 1 should be of double or int datatype;\
			arguments 2 and 3 should be of double datatype.\r\n", fname);
			return -1;
		}
		
		CvMat* _train_labels= SciMat2CvMat(1);
		CvMat* _train_data = SciMat2CvMat(2);
		CvMat* _test_data = SciMat2CvMat(3);
		CvMat* _test_labels = cvCreateMat(_test_data->rows, 1, CV_32FC1);
		
		int m, n, l;
		int _k;
		int _is_regression = 0;
		
		GetRhsVar(4, MATRIX_OF_DOUBLE_DATATYPE, &m, &n, &l);
		_k = (int)(*stk(l));
		
		if(Rhs == 5)
		{
			GetRhsVar(5, MATRIX_OF_DOUBLE_DATATYPE, &m, &n, &l);
			_is_regression = (int)*stk(l);
		}
		
		/////////////////////////////////////////////////////////////////////////////////////
	    //CvKNearest( const CvMat* _train_data, const CvMat* _responses,
	    //            const CvMat* _sample_idx=0, bool _is_regression=false, int max_k=32 );
		//
	    //virtual bool train( const CvMat* _train_data, const CvMat* _responses,
	    //                    const CvMat* _sample_idx=0, bool is_regression=false,
	    //                    int _max_k=32, bool _update_base=false );
		//
	    //virtual float find_nearest( const CvMat* _samples, int k, CvMat* results=0,
	    //    const float** neighbors=0, CvMat* neighbor_responses=0, CvMat* dist=0 ) const;
		/////////////////////////////////////////////////////////////////////////////////////
		
		CvMat* _sample_idx = 0;

		//learn classifier
		CvKNearest knn(_train_data, _train_labels, _sample_idx, _is_regression, _k);

		//float** neighbors = 0;
		//CvMat* neighbor_responses = 0;
		//CvMat* dist = 0;
		float response;
		
		response = knn.find_nearest(_test_data, _k, _test_labels, 0, 0, 0);

		//column-wise to row-wise
		//
		CvMat2SciMat(_test_labels, Rhs+1);
		
		LhsVar(1) = Rhs+1; 
		
		cvReleaseMat(&_train_data);
		cvReleaseMat(&_train_labels);
		cvReleaseMat(&_test_data);
		cvReleaseMat(&_test_labels);
		
		return 0;
	}
	

}


