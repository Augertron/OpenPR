///////////////////////////////////////////////////////////////////////////////
// Author:  Jia Wu
// Version: 0.1
// Date: Nov 2009
// Description: Generalized Discriminant Analysis(GDA)
// Reference: G. Baudat, F. Anouar, ¡°Generalized Discriminant Analysis Using
//            a Kernel Approach", Neural Computation, 12:2385-2404, 2000.
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
//        x          - dim*num data matrix. Each column is a data point.
//        c 	        - Class label vector of size 1*num or num*1.
//        ker        - A struct variable for creating kernel matrix.
//
//  Output:
//        norAlpha   - Alpha normalized vector.
///////////////////////////////////////////////////////////////////////////////

function norAlpha = generalizedda(x, c, ker)
  
  [dim, num] = size(x);
  [row, col] = size(c);
  
  tmp = (row+col-1)-num;
  if tmp,
    error('class label vector should has equal number of elements as sample number.');
  end
  
  if ~isstruct(ker),
    error('The third argument should be a struct variable.');
  end

  class_label = unique(c);
  class_num = length(class_label);
  
  [c, idx] = sort(c);
  x = x(:, idx);
  
  K = createkernel(x, x, ker);
  //center the kernel matrix
  K_n = K - ones(num, 1)*sum(K, 'r')/num - (ones(num, 1)*sum(K, 'r')/num)' + (sum(K)/num^2)*ones(num, num);
  
  [eigvec_K, eigval_K] = spec(K_n);
  
  [eigval_K, idx] = sort(diag(eigval_K));
  eigvec_K = eigvec_K(:, idx);
  
  //remove eigenvectors with very small eigenvalues
  thre = max(eigval_K)/1e4;
  idx = find(eigval_K > thre);
  eigvec_K = eigvec_K(:, idx);
  eigval_K = eigval_K(idx);
  rank_K = length(idx);
  
  //centered kernel matrix with only the selected large eigenvalues/eigenvectors
  K_n = eigvec_K * eigval_K * eigvec_K'; 
  
  //construct the block diagonal matrix W
  W = zeros(num, num);
  st = 1;
  ed = 0
  for i = 1:class_num,
    tmp = find(c==class_label(i));
    num_one_class = length(tmp);
    ed = ed+num_one_class;
    for j = st:ed,
      for k = st:ed,
        W(j, k) = 1/num_one_class;
      end
    end
    st = st+num_one_class;
  end

  new_dim = min((class_num-1), rank_K);

  //calculate normalized alpha 
  [eigvec_beta, eigval_beta] = spec(eigvec_K'*W*eigvec_K);
  [eigval_beta, idx] = sort(diag(eigval_beta));
  eigvec_beta = eigvec_beta(:, idx);
  eigvec_beta = eigvec_beta(:, 1:new_dim);
  
  alpha = eigvec_K*inv(eigval_K)*eigvec_beta;
  norAlpha = zeros(alpha);
  for i = 1:new_dim,
    norAlpha(:,i) = alpha(:,i)/sqrtm(alpha(:,i)'*K_n*alpha(:, i));   
  end
    
endfunction
