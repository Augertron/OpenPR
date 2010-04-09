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

#include "transformation.h"


#define NUM_OF_PARAMS_FIELDS 9
#define NUM_OF_MODEL_FIELDS 10
#define NDIM 3

const char *params_field_names[] = {
	"nclusters",
	"cov_mat_type",
	"start_step",
	"iter",
	"eps"
	"probs",
	"weights",
	"means",
	"covs"
};


const char *model_field_names[] = {
	"nclusters",
	"cov_mat_type",
	"log_likelihood",
	"means",
	"covs",
	"weights",
	"probs",
	"log_weights_div_det",
	"inv_eigen_values",
	"cov_rotate_mats"
};

void params_sci2cv(const mxArray *pArray, CvEMParams *params)
{
	//printf("params_sci2cv\n");

	int i, j, k, t;
	int nrow, ncol, ndim;

	int num_of_fields;
	int *dims;
	mxArray **rhs;
	int idx;
	double *ptr;

	//access "params"
	num_of_fields = mxGetNumberOfFields(pArray);
	if(num_of_fields != NUM_OF_PARAMS_FIELDS)
	{
		mexErrMsgTxt("params_sci2cv: number of returned 'params' fields is not correct.\n");
	}
	
	rhs = (mxArray **)mxMalloc(sizeof(mxArray *)*num_of_fields);
	for(i = 0; i < num_of_fields; i++)
	{
		rhs[i] = mxGetFieldByNumber(pArray, 0, i);
		//printf("params_sci2cv: get %dth field\n", i);
	}
	
	//nclusters
	idx = 0;
	ptr = mxGetPr(rhs[idx]);
	params->nclusters = (int)ptr[0];
	//printf("params_sci2cv - nclusters: %d\n", params->nclusters);
	idx++;
	
	//cov_mat_type
	ptr = mxGetPr(rhs[idx]);
	params->cov_mat_type = (int)ptr[0];
	//printf("params_sci2cv - cov_mat_type: %d\n", params->cov_mat_type);
	idx++;
	
	//start_step
	ptr = mxGetPr(rhs[idx]);
	params->start_step = (int)ptr[0];
	idx++;
	
/*	//term_crit 
	num_of_fields = mxGetNumberOfFields(rhs[idx]);
	if(num_of_fields != 2)
	{
		mexErrMsgTxt("params_sci2cv: term_crit expects 2 fields.\n");
	}	
	//iter
	mxArray *iter = mxGetFieldByNumber(rhs[idx], 0, 0);
	//eps
	mxArray *eps = mxGetFieldByNumber(rhs[idx], 0, 1);	
    params->term_crit = cvTermCriteria(CV_TERMCRIT_ITER+CV_TERMCRIT_EPS, mxGetPr(iter)[0], mxGetPr(eps)[0]);  
    idx++;*/
    
    //iter
   	ptr = mxGetPr(rhs[idx]);
	int max_iter = (int)ptr[0];
	//printf("params_sci2cv - iter\n");
	idx++;
	
	//eps
	ptr = mxGetPr(rhs[idx]);
	double epsilon = ptr[0];
	//printf("params_sci2cv - eps: %f\n", epsilon);	
	idx++;	
	
	params->term_crit = cvTermCriteria(CV_TERMCRIT_ITER+CV_TERMCRIT_EPS, max_iter, epsilon);
	
	//probs	
	if(!mxIsNumeric(rhs[idx]))
	{
		mexErrMsgTxt("params_sci2cv: probs must be numeric.\n");
	}	
	
	if(!mxIsEmpty(rhs[idx]))
	{
		nrow = mxGetM(rhs[idx]);
		ncol = mxGetN(rhs[idx]);	
		
		params->probs = cvCreateMat(nrow, ncol, CV_64FC1);
		
		ptr = mxGetPr(rhs[idx]);
		
		k = 0;
		for(j = 0; j < ncol; j++)
		{
			for(i = 0; i < nrow; i++)
			{
				((double*)(params->probs->data.ptr+params->probs->step*i))[j] = ptr[k++];
			}
		}
	}
	else
	{
		params->probs = NULL;
	}
	idx++;
	
	//weights
	if(!mxIsNumeric(rhs[idx]))
	{
		mexErrMsgTxt("params_sci2cv: weights must be numeric.\n");
	}
	
	if(!mxIsEmpty(rhs[idx]))
	{	
		nrow = mxGetM(rhs[idx]);
		ncol = mxGetN(rhs[idx]);
		
		params->weights = cvCreateMat(nrow, ncol, CV_64FC1);
		
		ptr = mxGetPr(rhs[idx]);
		
		ncol = ncol+nrow-1;		
		for(j = 0; j < ncol; j++)
		{
			params->weights->data.db[j] = ptr[j];
		}
	}
	else
	{
		params->weights = NULL;
	}
	idx++;
	
	//means
	if(!mxIsNumeric(rhs[idx]))
	{
		mexErrMsgTxt("params_sci2cv: means must be numeric.\n");
	}
	
	if(!mxIsEmpty(rhs[idx]))
	{
		nrow = mxGetM(rhs[idx]);
		ncol = mxGetN(rhs[idx]);
		
		params->means = cvCreateMat(nrow, ncol, CV_64FC1);
		
		ptr = mxGetPr(rhs[idx]);
		
		k = 0;
		for(j = 0; j < ncol; j++)
		{
			for(i = 0; i < nrow; i++)
			{
				((double*)(params->means->data.ptr+params->means->step*i))[j] = ptr[k++];
			}
		}
	}
	else
	{
		params->means = NULL;
	}
	idx++;
	
	//covs
	ptr = mxGetPr(rhs[idx]);
	
	if(!mxIsNumeric(rhs[idx]))
	{
		mexErrMsgTxt("params_sci2cv: covs must be numeric.\n");
	}	
	
	if(!mxIsEmpty(rhs[idx]))
	{
		ndim = mxGetNumberOfDimensions(rhs[idx]);
		if(ndim != NDIM)
		{
			mexErrMsgTxt("params_sci2cv: covs should be a 3-dimensional array.\n");	//CvMat** covs
		}
		
		dims = mxGetDimensions(rhs[idx]);
		nrow = dims[0];
		ncol = dims[1];
		
		params->covs = (const CvMat**)cvAlloc(dims[ndim]*sizeof(CvMat*));
				
		ptr = mxGetPr(rhs[idx]);

		k = 0;
		for(t = 0; t < dims[ndim]; t++)
		{
			params->covs[t] = cvCreateMat(nrow, ncol, CV_64FC1);
		
			for(j = 0; j < ncol; j++)
			{
				for(i = 0; i < nrow; i++)
				{
					((double*)(params->covs[t]->data.ptr+params->covs[t]->step*i))[j] = ptr[k++];
				}
			}		
		}
	}
	else
	{
		params->covs = NULL;
	}
    
    //release memory
    mxFree(rhs);
}

/*void params_cv2sci(CvEMParams *params, mxArray **pArray)
{
	//printf("params_cv2sci transfer\n");

	mxArray **lhs;
	
	int i, j, k, t;
	int idx;
	double *ptr;
	
	int nfields;
	int nrow, ncol;
	
	lhs = (mxArray**)mxMalloc(NUM_OF_PARAMS_FIELDS*sizeof(mxArray*));	
	//fill "params"
	//nclusters
	idx = 0;
	lhs[idx] = mxCreateDoubleMatrix(1, 1, mxREAL);
	ptr = mxGetPr(lhs[idx]);
	ptr[0] = (double)params->nclusters;
	idx++;	
	//printf("params-cv2sci transfer - nclusters\n");
	//printf("nclusters %f\n", ptr[0]);
	
	//cov_mat_type
	lhs[idx] = mxCreateDoubleMatrix(1, 1, mxREAL);
	ptr = mxGetPr(lhs[idx]);
	ptr[0] = (double)params->cov_mat_type;
	idx++;
	//printf("params-cv2sci transfer - cov_mat_type\n");
	//printf("cov_mat_type %f\n", ptr[0]);	
	
	//start_step
	lhs[idx] = mxCreateDoubleMatrix(1, 1, mxREAL);
	ptr = mxGetPr(lhs[idx]);
	ptr[0] = (double)params->start_step;
	idx++;
	//printf("params-cv2sci transfer - start_step\n");
	//printf("start_step %f\n", ptr[0]);	
	
	//term_crit
	nfields = 2;
	const char *field_names[] = {"iter", "eps"};
	
	lhs[idx] = mxCreateStructMatrix(1, 1, nfields, field_names);
	double val = params->term_crit.max_iter;
	mxSetFieldByNumber(lhs[idx], 0, 0, mxCreateDoubleScalar(val));
	val = params->term_crit.epsilon;
	mxSetFieldByNumber(lhs[idx], 0, 1, mxCreateDoubleScalar(val));	
	idx++;
	//printf("params-cv2sci transfer - term_crit\n");
	
	//probs
	if(!params->probs)
	{
		lhs[idx] = mxCreateDoubleMatrix(0, 0, mxREAL);
	}
	else
	{
		nrow = params->probs->rows;
		ncol = params->probs->cols;
		
		lhs[idx] = mxCreateDoubleMatrix(nrow, ncol, mxREAL);		
		ptr = mxGetPr(lhs[idx]);
		
		k = 0;
		for(j = 0; j < ncol; j++)
		{
			for(i = 0; i < nrow; i++)
			{
				ptr[k++] = CV_MAT_ELEM(*params->probs, double, i, j);
			}
		}
	}
	idx++;
	//printf("params-cv2sci transfer - probs\n");
	
	//weights
	if(!params->weights)
	{
		lhs[idx] = mxCreateDoubleMatrix(0, 0, mxREAL);
	}
	else
	{
		nrow = params->weights->rows;
		ncol = params->weights->cols;
		
		lhs[idx] = mxCreateDoubleMatrix(nrow, ncol, mxREAL);
		ptr = mxGetPr(lhs[idx]);
		
		k = 0;
		for(j = 0; j < ncol; j++)
		{
			for(i = 0; i < nrow; i++)
			{
				ptr[k++] = CV_MAT_ELEM(*params->weights, double, i, j);
			}
		}		
	}
	idx++;
	//printf("params-cv2sci transfer - weights\n");
	
	//means
	if(!params->means)
	{
		lhs[idx] = mxCreateDoubleMatrix(0, 0, mxREAL);
	}
	else
	{
		nrow = params->means->rows;
		ncol = params->means->cols;
		
		lhs[idx] = mxCreateDoubleMatrix(nrow, ncol, mxREAL);
		ptr = mxGetPr(lhs[idx]);
		
		k = 0;
		for(j = 0; j < ncol; j++)
		{
			for(i = 0; i < nrow; i++)
			{
				ptr[k++] = CV_MAT_ELEM(*params->means, double, i, j);
			}
		}
	}
	idx++;
	//printf("params-cv2sci transfer - means\n");
	
	//covs
	if(!params->covs)
	{	
		lhs[idx] = mxCreateDoubleMatrix(0, 0, mxREAL);		
	}
	else
	{
		const int dims[]={params->covs[0]->rows, params->covs[0]->cols, params->nclusters};
		
		lhs[idx] = mxCreateNumericArray(NDIM, dims, mxDOUBLE_CLASS, mxREAL);
		ptr = mxGetPr(lhs[idx]);
		
		k = 0;
		for(t = 0; t < dims[2]; t++)
		{
			for(j = 0; j < dims[1]; j++)
			{
				for( i = 0; i < dims[0]; i++)
				{
					ptr[k++] = CV_MAT_ELEM(*params->covs[t], double, i, j);
				}
			}
		}
	}	
	//printf("params-cv2sci transfer - covs\n");
	
	//return "params"
	pArray[0] = mxCreateStructMatrix(1, 1, NUM_OF_PARAMS_FIELDS, params_field_names);
	for(i = 0; i < NUM_OF_PARAMS_FIELDS; i++)
	{
		mxSetField(pArray[0], 0, params_field_names[i], lhs[i]);
		//printf("set field %d name %s\n", i, params_field_names[i]);
	}
	
	//release memory
	mxFree(lhs);	
	//printf("params-cv2sci free mem\n");
}*/

void model_sci2cv(const mxArray *pArray, CvEM *em_model)
{
	//printf("model_sci2cv\n");	
	
	mxArray **rhs;
	
	int num_of_fields;
	int idx;
	double *ptr;
	
	int nrow, ncol;
	int ndim;
	int *dims = NULL;
	int i, j, k, t;
	
	if(!em_model)
	{
		mexErrMsgTxt("model_sci2cv: object of CvEM does not exist.\n");
	}
	
	//access "model"
	num_of_fields = mxGetNumberOfFields(pArray);
	if(num_of_fields != NUM_OF_MODEL_FIELDS)
	{
		mexErrMsgTxt("model_sci2cv: number of returned 'model' fields is not correct.\n");
	}
	
	rhs = (mxArray**)mxMalloc(NUM_OF_MODEL_FIELDS*sizeof(mxArray*));
	
	for(i = 0; i < num_of_fields; i++)
	{
		rhs[i] = mxGetFieldByNumber(pArray, 0, i);
	}
	
/*	//params
	idx = 0;
	params_sci2cv(rhs[idx], &(em_model->params));
	idx++;
	//printf("model_sci2cv: access 'params'\n");*/
	
	//fill CvEM object
	//nclusters
	idx = 0;
	ptr = mxGetPr(rhs[idx]);
	em_model->params.nclusters = (int)ptr[0];
	idx++;
	//printf("model_sci2cv - nclusters\n");	
	
	//cov_mat_type
	ptr = mxGetPr(rhs[idx]);
	em_model->params.cov_mat_type = (int)ptr[0];
	idx++;
	//printf("model_sci2cv - cov_mat_type\n");
	
	//log_likelihood
	ptr = mxGetPr(rhs[idx]);
	em_model->log_likelihood = ptr[0];
	idx++;
	//printf("model_sci2cv - log_likelihood\n");	
	
	//means
	nrow = mxGetM(rhs[idx]);
	ncol = mxGetN(rhs[idx]);
	em_model->means = cvCreateMat(nrow, ncol, CV_64FC1);
	ptr = mxGetPr(rhs[idx]);
	k = 0;	
	for(j = 0; j < ncol; j++)
	{
		for(i = 0; i < nrow; i++)
		{
			CV_MAT_ELEM(*(em_model->means), double, i, j) = ptr[k++];
		}
	}
	idx++;
	//printf("model_sci2cv - means\n");	
	
	//covs
	dims = mxGetDimensions(rhs[idx]);
	nrow = dims[0];
	ncol = dims[1];		
	ndim = dims[2];
	
	em_model->covs = (CvMat**)cvAlloc(ndim*sizeof(CvMat*));
	ptr = mxGetPr(rhs[idx]);
	
	k = 0;
	for(t = 0; t < ndim; t++)
	{
		em_model->covs[t] = cvCreateMat(nrow, ncol, CV_64FC1);
		
		for(j = 0; j < ncol; j++)
		{
			for(i = 0; i < nrow; i++)
			{
				CV_MAT_ELEM(*(em_model->covs[t]), double, i, j) = ptr[k++];
			}
		}
	}
	idx++;
	//printf("model_sci2cv - covs\n");
	
	//weights
	nrow = mxGetM(rhs[idx]);
	ncol = mxGetN(rhs[idx]);
	em_model->weights = cvCreateMat(nrow, ncol, CV_64FC1);
	ptr = mxGetPr(rhs[idx]);
	k = 0;
	for(j = 0; j < ncol; j++)
	{
		for(i = 0; i < nrow; i++)
		{
			CV_MAT_ELEM(*(em_model->weights), double, i, j) = ptr[k++];
		}
	}
	idx++;
	//printf("model_sci2cv - weights\n");	
	
	//probs
	nrow = mxGetM(rhs[idx]);
	ncol = mxGetN(rhs[idx]);
	em_model->probs = cvCreateMat(nrow, ncol, CV_64FC1);
	ptr = mxGetPr(rhs[idx]);
	k = 0;
	for(j = 0; j < ncol; j++)
	{
		for(i = 0; i < nrow; i++)
		{
			CV_MAT_ELEM(*(em_model->probs), double, i, j) = ptr[k++];
		}
	}
	idx++;
	//printf("model_sci2cv - probs\n");
	
	//log_weight_div_det
	nrow = mxGetM(rhs[idx]);
	ncol = mxGetN(rhs[idx]);
	em_model->log_weight_div_det = cvCreateMat(nrow, ncol, CV_64FC1);
	ptr = mxGetPr(rhs[idx]);
	k = 0;
	for(j = 0; j < ncol; j++)
	{
		for(i = 0; i < nrow; i++)
		{
			CV_MAT_ELEM(*(em_model->log_weight_div_det), double, i, j) = ptr[k++];
		}
	}
	idx++;
	//printf("model_sci2cv - log_weight_div_det\n");	
	
	//inv_eigen_values
	nrow = mxGetM(rhs[idx]);
	ncol = mxGetN(rhs[idx]);
	em_model->inv_eigen_values = cvCreateMat(nrow, ncol, CV_64FC1);
	ptr = mxGetPr(rhs[idx]);
	k = 0;
	for(j = 0; j < ncol; j++)
	{
		for(i = 0; i < nrow; i++)
		{
			CV_MAT_ELEM(*(em_model->inv_eigen_values), double, i, j) = ptr[k++];
		}
	}
	idx++;
	//printf("model_sci2cv - inv_eigen_values\n");
		
	//cov_rotate_mats
	dims = mxGetDimensions(rhs[idx]);
	nrow = dims[0];
	ncol = dims[1];	
	ndim = dims[2];
	
	em_model->cov_rotate_mats = (CvMat**)cvAlloc(ndim*sizeof(CvMat*));
	ptr = mxGetPr(rhs[idx]);
	
	k = 0;
	for(t = 0; t < ndim; t++)
	{
		em_model->cov_rotate_mats[t] = cvCreateMat(nrow, ncol, CV_64FC1);
		for(j = 0; j < ncol; j++)
		{
			for(i = 0; i < nrow; i++)
			{
				CV_MAT_ELEM(*(em_model->cov_rotate_mats[t]), double, i, j) = ptr[k++];
			}
		}
	}
	//printf("model_sci2cv - cov_rotate_mats\n");
	
	//release memory
	mxFree(rhs);
}

void model_cv2sci(CvEM *em_model, mxArray **pArray)
{
	//printf("model_cv2sci transfer\n");	
	
	mxArray **rhs = (mxArray**)mxMalloc(NUM_OF_MODEL_FIELDS*sizeof(mxArray*));
	mxArray *return_model;
	
	int idx;
	double *ptr;
	int nrow, ncol;
	int i, j, k, t;
	
/*	//params
	idx = 0;
	params_cv2sci(&(em_model->params), &rhs[idx]);
	idx++;
	//printf("model_cv2sci transfer - params\n");*/
	
	//nclusters
	idx = 0;
	rhs[idx] = mxCreateDoubleMatrix(1, 1, mxREAL);
	ptr = mxGetPr(rhs[idx]);
	ptr[0] = (double)em_model->params.nclusters;
	idx++;
	
	//cov_mat_type
	rhs[idx] = mxCreateDoubleMatrix(1, 1, mxREAL);
	ptr = mxGetPr(rhs[idx]);
	ptr[0] = (double)em_model->params.cov_mat_type;
	//printf("model_cv2sci transfer - cov_mat_type\n");
	//printf("cov_mat_type: %d\n", em_model->params.cov_mat_type);
	idx++;
	
	//log_likelihood
	rhs[idx] = mxCreateDoubleMatrix(1, 1, mxREAL);
	ptr = mxGetPr(rhs[idx]);
	ptr[0] = em_model->log_likelihood;
	idx++;
	//printf("model_cv2sci transfer - log_likelihood\n");
	
	//means
	nrow = em_model->means->rows;
	ncol = em_model->means->cols;
	rhs[idx] = mxCreateDoubleMatrix(nrow, ncol, mxREAL);
	ptr = mxGetPr(rhs[idx]);
	
	k = 0;
	for(j = 0; j < ncol; j++)
	{
		for(i = 0; i < nrow; i++)
		{
			ptr[k++] = CV_MAT_ELEM(*em_model->means, double, i, j);
		}
	}
	idx++;
	//printf("model_cv2sci transfer - means\n");
	
	//covs
	int dims[]={em_model->covs[0]->rows, em_model->covs[0]->cols, em_model->params.nclusters};	
	
	rhs[idx] = mxCreateNumericArray(NDIM, dims, mxDOUBLE_CLASS, mxREAL);
	ptr = mxGetPr(rhs[idx]);
	
	k = 0;
	for(t = 0; t < dims[2]; t++)
	{
		for(j = 0; j < dims[1]; j++)
		{
			for(i = 0; i < dims[0]; i++)
			{
				ptr[k++] = CV_MAT_ELEM(*em_model->covs[t], double, i, j);
			}
		}
	}
	idx++;
	//printf("model_cv2sci transfer - covs\n");
	
	//weights
	nrow = em_model->weights->rows;
	ncol = em_model->weights->cols;
	
	rhs[idx] = mxCreateDoubleMatrix(nrow, ncol, mxREAL);
	ptr = mxGetPr(rhs[idx]);
	
	k = 0;
	for(j = 0; j < ncol; j++)
	{
		for(i = 0; i < nrow; i++)
		{
			ptr[k++] = CV_MAT_ELEM(*em_model->weights, double, i, j);
		}
	}
	idx++;
	//printf("model_cv2sci transfer - weights\n");
	
	//probs
	nrow = em_model->probs->rows;
	ncol = em_model->probs->cols;
	
	rhs[idx] = mxCreateDoubleMatrix(nrow, ncol, mxREAL);
	ptr = mxGetPr(rhs[idx]);
	
	k = 0;
	for(j = 0; j < ncol; j++)
	{
		for(i = 0; i < nrow; i++)
		{
			ptr[k++] = CV_MAT_ELEM(*em_model->probs, double, i, j);
			//printf("%.10f ", CV_MAT_ELEM(*em_model->probs, double, i, j));
		}
		//printf("\n");
	}
	idx++;
	//printf("model_cv2sci transfer - probs\n");
	
	//log_weight_div_det
	nrow = em_model->log_weight_div_det->rows;
	ncol = em_model->log_weight_div_det->cols;
	
	rhs[idx] = mxCreateDoubleMatrix(nrow, ncol, mxREAL);
	ptr = mxGetPr(rhs[idx]);
	
	k = 0;
	for(j = 0; j < ncol; j++)
	{
		for(i = 0; i < nrow; i++)
		{
			ptr[k++] = CV_MAT_ELEM(*em_model->log_weight_div_det, double, i, j);
		}
	}
	idx++;
	//printf("model_cv2sci transfer - log_weight_div_det\n");
	
	//inv_eigen_values
	nrow = em_model->inv_eigen_values->rows;
	ncol = em_model->inv_eigen_values->cols;
	
	rhs[idx] = mxCreateDoubleMatrix(nrow, ncol, mxREAL);
	ptr = mxGetPr(rhs[idx]);
	
	k = 0;
	for(j = 0; j < ncol; j++)
	{
		for(i = 0; i < nrow; i++)
		{
			ptr[k++] = CV_MAT_ELEM(*em_model->inv_eigen_values, double, i, j);
		}
	}
	idx++;
	//printf("model_cv2sci transfer - inv_eigen_values\n");
	
	//cov_rotate_mats
	dims[0] = em_model->cov_rotate_mats[0]->rows;
	dims[1] = em_model->cov_rotate_mats[0]->cols;
	dims[2] = em_model->params.nclusters;
	
	rhs[idx] = mxCreateNumericArray(NDIM, dims, mxDOUBLE_CLASS, mxREAL);
	ptr = mxGetPr(rhs[idx]);
			
	k = 0;
	for(t = 0; t < dims[2]; t++)
	{
		for(j = 0; j < dims[1]; j++)
		{
			for(i = 0; i < dims[0]; i++)
			{
				ptr[k++] = CV_MAT_ELEM(*em_model->cov_rotate_mats[t], double, i, j);
			}
		}
	}
	//printf("model_cv2sci transfer - cov_rotate_mats\n");	
	
	//write and return model
	return_model = mxCreateStructMatrix(1, 1, NUM_OF_MODEL_FIELDS, model_field_names);
	for(i = 0; i < NUM_OF_MODEL_FIELDS; i++)
	{
		mxSetField(return_model, 0, model_field_names[i], mxDuplicateArray(rhs[i]));
		//printf("set field %d name %s\n", i, model_field_names[i]);
	}
	pArray[0] = return_model;
	
	//release memory
	mxFree(rhs);
}
