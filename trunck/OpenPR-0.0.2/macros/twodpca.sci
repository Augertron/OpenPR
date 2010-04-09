///////////////////////////////////////////////////////////////////////////////
// Author:  Jia Wu
// Version: 0.1
// Date: Nov 2009
// Description: Two-Dimensional PCA
// Reference:
//          J. Yang, D. Zhang, A. Frangi, and J. Yang. Two-dimensional pca: A 
//          new approach to appearance-based face representation and recognition. 
//          IEEE Trans. on Pattern Analysis and Machine Intelligence, 
//          26(1):131–137, 2004.
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
/////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//  Input:
//        A       - m*n*d hypermat
//				    m*n is the size of the sample matrix
//                  d is the number of the sample
//		  K       - number of projection axis Xk; K <= n
//
//  Output:
//        Y       - m*K*d hypermat
//                  m*K is the size of the feature matrix; each column is
//                  a projected feature vector Yk = A*Xk
//				    d is the number of the sample
//		  X		  - The set of projection axes; each column is a projection axis.
///////////////////////////////////////////////////////////////////////////////

function [Y, X] = twodpca(A, K)
	
	[m, n, d] = size(A);
	
	if(argn(2) < 2),
		K = n;
	end
	if(K < 0)
		error('K should be positive.');
	end
	if(K > n)
		error('K should be no more than number of columns of the sample matrix');
	end
	
	Aavg = zeros(m, n);
	for i = 1:d,
		Aavg = Aavg+A(:,:,i);
	end
	Aavg = Aavg./d;
	
	Gt = zeros(n, n);
	for j = 1:d,
		Gt = Gt+((A(:,:,i)-Aavg)')*(A(:,:,i)-Aavg);
	end
	Gt = Gt./d;
	
	[eig_vec, eig_val] = spec(Gt);
	[eig_val, idx] = sort(diag(eig_val));
	eig_vec = eig_vec(:, idx);
	
	X = eig_vec(:, 1:K);
	
	for t = 1:d,
		Y(:,:,t) = A(:,:,t)*X;
	end

endfunction
