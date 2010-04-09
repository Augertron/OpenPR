///////////////////////////////////////////////////////////////////////////////
// Author:  Jia Wu
// Version: 0.1
// Date: Nov 2009
// Description: Kernel Principal Component Analysis
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
//        ker            - The kernel matrix.
//		  dimension		 - Number of dimension for the new patterns.
//
//  Output:
//        new_patterns   - New patterns.
//        eig_val        - The sorted eigenvalue.
//        eig_vec        - Sorted eigenvector
//                         Each column is an eigenvector. eig_vec'*eig_vec=I.
///////////////////////////////////////////////////////////////////////////////

function [eig_vec, eig_val, new_patterns] = kpca(patterns, ker, dimension)
  
  K = createkernel(patterns, [], ker);
  
  [r, c] = size(K);
  num = size(patterns, 2);  
  
  ln = 1/num*ones(r, c); 
  Kn = K-ln*K-K*ln+ln*K*ln;
  
  [eig_vec, eig_val] = spec(Kn);
  [eig_val, idx] = sort(diag(eig_val));
  eig_vec = eig_vec(:, idx);
  
  eig_vec = eig_vec(:, 1:dimension);
  
  new_patterns = eig_vec'*K;  
  
endfunction

