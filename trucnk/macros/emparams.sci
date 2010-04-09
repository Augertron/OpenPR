//////////////////////////////////////////////////////////////////////////////
// Author:  Jia Wu
// Version: 0.1
// Date: Dec. 2009
//
//
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

function params = emparams(nclusters, cov_mat_type, start_step, iter, eps, probs, weights, means, covs)
//function params = emparams(nclusters, cov_mat_type, start_step, term_crit, probs, weights, means, covs)

	rhs = argn(2);
	
	if rhs == 0 then
		disp('Create a struct variable containing parameter values for EM algorithm.');
		disp('Usage:', 'params = emparams(nclusters[,cov_mat_type,start_step,iter,eps,probs,weights,means,covs])');
		disp('nclusters: number of Gaussian distributions.');
		disp('cov_mat_type: type of covariance matrix. 0 - spherical, 1 - diagonal, 2 - generic');
		disp('start_step: initial step EM starts from. 0 - Auto-step, 1 - E-step, 2 - M-step');
		disp('iter/eps: termination criteria of the procedure. iter for iteration times; eps for difference of change.');			
		disp('probs: initial probabilities (Pi,k); used(must be not NULL) only when EM starts from M-step');
		disp('weights: initial weights for each distribution; used(if not NULL) only when EM starts from E-step');
		disp('means: initial means of each distribution; used(must be not NULL) only when EM starts from E-step');
		disp('covs: initial covariance matrix of each distribution; used(if not NULL) only when EM starts from E-step');
//		disp('term_crit: termination criteria of the procedure. It is struct variable and is a combination of iteration times and parameter of change');	
		disp('params: a struct variable containing the parameters of EM algorithm.');	
		abort;
	end
	
	if rhs < 1 | rhs > 8 then
		error('The number of input arguments should be in the range of [1, 8].\n');
	end
	
	//default
	cov_mat_type1 = 1;
	start_step1 = 0;
	probs1 = [];
	weights1 = [];
	means1 = [];
	covs1 = [];
//	term_crit1 = struct('iter', 1000, 'eps', 1e-6);
	iter1 = 1000;
	eps1 = 1e-6;
	
	if exists('cov_mat_type') then
		cov_mat_type1 = cov_mat_type;
	end
	if exists('start_step') then
		start_step1 = start_step; 
	end
//	if exists('term_crit') then
//		if ~isstruct(term_crit) then
//			error('The argument term_crit must be a struct variable.');
//		end
		
//		term_crit1 = term_crit;
//	end
	if exists('iter') then
		iter1 = iter;
	end
	if exists('eps') then
		eps1 = eps;
	end
	if exists('probs') then
		probs1 = probs;
	end
	if exists('weights') then
		weights1 = weights;
	end 
	if exists('means') then
		means1 = means;
	end
	if exists('covs') then
		ndim = length(size(covs));
		if ndim ~= 3 then
			error('The argument covs should be a 3-dimensional array.');
		end
		
		nc = size(covs, 3);
		if nc ~= nclusters then
			error('The sizes of nclusters and covs do not match.');
		end
		
		covs1 = covs;
	end
			
//	params = struct('nclusters',nclusters, 'cov_mat_type',cov_mat_type1, 'start_step',start_step1, 'term_crit',term_crit1, 'probs',probs1, 'weights',weights1, 'means',means1, 'covs',covs1);

	params = struct('nclusters',nclusters, 'cov_mat_type',cov_mat_type1, 'start_step',start_step1, 'iter',iter1, 'eps',eps1, 'probs',probs1, 'weights',weights1, 'means',means1, 'covs',covs1);

endfunction
