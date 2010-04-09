///////////////////////////////////////////////////////////////////////////////
// Author:  Jia Wu
// Version: 0.1
// Date: Nov 2009
// Description: Principal Component Analysis(PCA)
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
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//  Input:
//        patterns       - Data matrix. Each column is a data point.
//		  dimension		 - Number of dimension for the new patterns.
//
//  Output:
//        new_patterns   - New patterns.
//		  m				 - Mean. 
//        eig_val        - The sorted eigenvalue.
//        eig_vec        - Each column is an eigenvector. eig_vec'*eig_vec=I.
///////////////////////////////////////////////////////////////////////////////


function [new_patterns, m, eig_val, eig_vec] = pca2(patterns, dimension)

	[dim, sample_num] = size(patterns);
	
	if(dim < dimension),
		disp('Final dimension must not be larger than the data dimension.')
		disp('Set dimension to be the data dimension.')
		dimension = dim;
	end
	
	//mean
	m = mean(patterns, 'c');
	//covariance matrix
	cov_mat = (1/(sample_num-1))*(patterns-m*ones(1,sample_num))*(patterns-m*ones(1,sample_num))';
	//eigenvalues and eigenvecotrs
	[eig_vec, eig_val] = spec(cov_mat);
	
	//sort eigenvectors
	[eig_val, I] = sort(diag(eig_val));
	eig_vec = eig_vec(:, I);
	
	V = eig_vec(:, 1:dimension);	
	//new patterns
	new_patterns = V'*(patterns-m*ones(1,sample_num));

endfunction
