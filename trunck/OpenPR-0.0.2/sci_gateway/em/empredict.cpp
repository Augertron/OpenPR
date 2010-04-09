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
	#define mexFunction mex_empredict
#endif


#include "cxcore.h"

#include "ml.h"

#include "mex.h"

#include "transformation.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if(nrhs == 0)
	{
		mexPrintf("Usage:\n\t[labels, probs] = empredict(test_samples, em_model).\n");
		return;
	}
	
	//check arguments
	if(nrhs != 2)
	{
		mexErrMsgTxt("Expects 2 input arguments.");
	}
	
	if(nlhs != 2)
	{
		mexErrMsgTxt("Expects 2 output arguments.");
	}	
	
	int nrow, ncol;
	int nclusters;
	int i, j, k;
	
	int ndim;
	double *ptr;
	int *pi;
	
	CvEM *em_model = new CvEM();
	CvMat **test_samples;
	
	CvMat *labels;
	CvMat **probs;
	
	//acess "test_samples"
	ndim = mxGetNumberOfDimensions(prhs[0]);
	if(ndim > 2)
	{
		mexErrMsgTxt("Error: test_samples must be a row vector or a matrix containing one sample per row.");
	}
	
	nrow = mxGetM(prhs[0]);
	ncol = mxGetN(prhs[0]);
	
	ptr = mxGetPr(prhs[0]);
	
	labels = cvCreateMat(1, nrow, CV_32SC1);
	
	test_samples = (CvMat**)cvAlloc(nrow*sizeof(CvMat*));
	for(i = 0; i < nrow; i++)
	{
		test_samples[i] = cvCreateMat(1, ncol, CV_32FC1);
	}
	k = 0;
	for(j = 0; j < ncol; j++)
	{
		for(i = 0; i < nrow; i++)
		{
			CV_MAT_ELEM(*(test_samples[i]), float, 0, j) = (float)ptr[k++];
		}
	}		
	//printf("empredict: access 'test_samples'\n");
	
	//acess "em_model"
	//printf("empredict: access 'em_model'\n");	
	model_sci2cv(prhs[1], em_model);
	nclusters = em_model->get_nclusters();
	
/*	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//check "em_model"
		
	printf("empredict: nclusters %d\n", nclusters);
	printf("empredict: cov_mat_type %d\n", em_model->params.cov_mat_type);
	printf("empredict: log_likelihood %f\n", em_model->log_likelihood);
	printf("empredict: means\n");
	const CvMat* pmeans = em_model->get_means();
	for(i = 0; i < pmeans->rows; i++)
	{
		for(j = 0; j < pmeans->cols; j++)
		{
			printf("%f ", CV_MAT_ELEM(*pmeans, double, i, j));
		}
		printf("\n");
	}
	printf("empredict: covs\n");
	const CvMat** pcovs = em_model->get_covs();
	for(k = 0; k < nclusters; k++)
	{
		for(i = 0; i < pcovs[0]->rows ; i++)
		{
			for(j = 0; j < pcovs[0]->cols; j++)
			{
				printf("%f ", CV_MAT_ELEM(*pcovs[k], double, i, j));
			}
			printf("\n");
		}
		printf("\n");
	}
	/*printf("empredict: weights\n");
	const CvMat* pweights = em_model->get_weights();
	for(i = 0; i < pweights->rows; i++)
	{
		for(j = 0; j < pweights->cols; j++)
		{
			printf("%f ", CV_MAT_ELEM(*pweights, double, i, j));
		}
		printf("\n");
	}*/
/*	printf("empredict: probs\n");
	const CvMat* pprobs = em_model->get_probs();
	for(i = 0; i < pprobs->rows; i++)
	{
		for(j = 0; j < pprobs->cols; j++)
		{
			printf("%f ", CV_MAT_ELEM(*pprobs, double, i, j));
		}
		printf("\n");
	}
	printf("empredict: log_weight_div_det\n");
	const CvMat* plog_weight_div_det = em_model->get_log_weight_div_det();
	for(i = 0; i < plog_weight_div_det->rows; i++)
	{
		for(j = 0; j < plog_weight_div_det->cols; j++)
		{
			printf("%f ", CV_MAT_ELEM(*plog_weight_div_det, double, i, j));
		}
		printf("\n");
	}	
	printf("empredict: inv_eigen_values\n");
	const CvMat* pinv_eigen_values = em_model->get_inv_eigen_values();
	for(i = 0; i < pinv_eigen_values->rows; i++)
	{
		for(j = 0; j < pinv_eigen_values->cols; j++)
		{
			printf("%f ", CV_MAT_ELEM(*pinv_eigen_values, double, i, j));
		}
		printf("\n");
	}		
	printf("empredict: cov_rotate_mats\n");
	const CvMat** pcov_rotate_mats = em_model->get_cov_rotate_mats();
	for(k = 0; k < nclusters; k++)
	{
		for(i = 0; i < pcov_rotate_mats[0]->rows ; i++)
		{
			for(j = 0; j < pcov_rotate_mats[0]->cols; j++)
			{
				printf("%f ", CV_MAT_ELEM(*pcov_rotate_mats[k], double, i, j));
			}
			printf("\n");
		}
		printf("\n");
	}	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	
	probs = (CvMat**)cvAlloc(nrow*sizeof(CvMat*));
	for(i = 0; i < nrow; i++)
	{
		probs[i] = cvCreateMat(1, nclusters, CV_64FC1);
	}
	
	//EM prediction
	int temp;
	for(i = 0; i < nrow; i++)
	{
		temp = cvRound(em_model->predict(test_samples[i], probs[i]));
		CV_MAT_ELEM(*labels, int, 0, i) = temp;

		//printf("samples[%d] %f %f\n", i, test_samples[i]->data.fl[0], test_samples[i]->data.fl[1]);
		//printf("labels%d: %d\n", i, temp);
	}
	//printf("empredict: prediction over\n");
	
	//write "labels"
	plhs[0] = mxCreateNumericMatrix(1, nrow, mxINT32_CLASS, mxREAL);
	pi = (int*)mxGetData(plhs[0]);
	
	k = 0;
	for(i = 0; i < nrow; i++)
	{
		pi[k] = CV_MAT_ELEM(*labels, int, 0, i);
		k++;
	}
	//printf("empredict: write 'labels'\n");
	
	//write "probs"
	plhs[1] = mxCreateDoubleMatrix(nrow, nclusters, mxREAL);
	ptr = mxGetPr(plhs[1]);
	
	k = 0;
	for(j = 0; j < nclusters; j++)
	{
		for(i = 0; i < nrow; i++)
		{
			ptr[k] = CV_MAT_ELEM(*probs[i], double, 0, j);
			//printf("%.10f ", ptr[k]);
			k++;
//			//printf("%.10f ", CV_MAT_ELEM(*probs[i], double, 0, j));
		}	
		//printf("\n");
	}	
	//printf("empredict: write 'probs'\n");
	
	//release memory
	cvReleaseMat(&labels);
	
	for(i = 0; i < nrow; i++)
	{
		cvReleaseMat(&test_samples[i]);
		cvReleaseMat(&probs[i]);
	}
	cvFree(&test_samples);
	cvFree(&probs);
	//printf("empredict: release memory\n");
}
