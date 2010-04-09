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

#ifdef _MSC_VER
	#define mexFunction mex_emtrain
#endif

#include "cxcore.h"
#include "ml.h"

#include "mex.h"

#include "transformation.h"

#define NUM_OF_PARAMS_FIELDS 9


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if(nrhs == 0)
	{
		mexPrintf("Usage:\n\t[model, labels] = emtrain(samples, params)\n");
		return;
	}
	
	if(nrhs != 2)
	{
		mexErrMsgTxt("Expects 2 input arguments.");
	}
	
	if(nlhs != 2)
	{
		mexErrMsgTxt("Expects 2 output arguments.");
	}
	
	CvEM em_model;
	CvMat *samples;
	CvMat *labels;
	CvMat *sample_idx = 0;
	CvEMParams params;	
	
	int i, j, k;
	
	double *samples_tmp;
	int nrow, ncol;	
	int num_of_fields;
	
	int *pi;
	
	//check arguments
	if(!mxIsDouble(prhs[0]))
	{
		mexErrMsgTxt("Input argument 1 must be a matrix containing samples one per row.\n");
	}
	
	if(!mxIsStruct(prhs[1]))
	{
		mexErrMsgTxt("Input argument 2 must be a struct variable.\n");
	}
	num_of_fields = mxGetNumberOfFields(prhs[1]);
	if(num_of_fields != NUM_OF_PARAMS_FIELDS)
	{
		mexErrMsgTxt("The number of fields of argument 2 does not match the needed value.\n");
	}
	
	//access "samples"
	samples_tmp = mxGetPr(prhs[0]);
	nrow = mxGetM(prhs[0]);
	ncol = mxGetN(prhs[0]);
	
	samples = cvCreateMat(nrow, ncol, CV_32FC1);
	k = 0;
	for(j = 0; j < ncol; j++)
	{
		for(i = 0; i < nrow; i++)
		{
			CV_MAT_ELEM(*samples, float, i, j) = (float)samples_tmp[k++];
		}
	}
	//printf("emtrain access - samples\n");
	
	//aceess "params"
	params_sci2cv(prhs[1], &params);
	//printf("emtrain access - params\n");
	//printf("emtrain nclusters: %d\n", params.nclusters);
	//printf("emtrain cov_mat_type: %d\n", params.cov_mat_type);
	//printf("emtrain eps: %f\n", params.term_crit.epsilon);
    
    //trainning
    labels = cvCreateMat(1, samples->rows, CV_32SC1);
    em_model.train(samples, sample_idx, params, labels);    
    /*const CvMat *pProbs = em_model.get_probs();
    for(i = 0; i < pProbs->rows; i++)
    {
    	for(j = 0; j < pProbs->cols; j++)
    	{
    		printf("%f ", CV_MAT_ELEM(*pProbs, double, i, j));
    	}
    	printf("\n");
    }*/
    //printf("training ends\n");
    
    //write "model"
    model_cv2sci(&em_model, &plhs[0]);
    //printf("emtrain write - model\n");
    
    //write "labels"
    nrow = labels->rows;
    ncol = labels->cols;
    
    plhs[1] = mxCreateNumericMatrix(nrow, ncol, mxINT32_CLASS, mxREAL);
    pi = (int*)mxGetData(plhs[1]);
    
    k = 0;
    for(j = 0; j < ncol; j++)
    {
    	for(i = 0; i < nrow; i++)
    	{
    		pi[k] = CV_MAT_ELEM(*labels, int, i, j);
   		    k++;
    	}
    }   
    //printf("emtrain write - labels\n"); 
    
    //free memory
    cvReleaseMat(&samples);
    cvReleaseMat(&labels);
}
