///////////////////////////////////////////////////////////////////////////////
// Author:  Jia Wu
// Version: 0.1
// Date: Nov 2009
// Description: Kernel Fisher Discriminant Analysis
// Reference: S. Mika, G. Ratsch, J. Weston, B. Scholkopf, and K.-R. M¨¹ller,
//            ¡°Fisher discriminant analysis with kernels,¡± in Neural Networks for
//            Signal Processing IX, Y.-H. Hu, J. Larsen, E. Wilson, and S. Douglas,
//            Eds. Piscataway, NJ: IEEE, 1999, pp. 41¨C48.
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
//  Note the function is used for two-class problems.
//  Input:
//        x       - Data matrix of size dim*num. Each column is a data point.
//    	  c       - Class label vector of size 1*num or num*1.
//        ker     - Kernel, a struct variable.
//		
//  Output:
//        alpha   - The principal eigenvector.
///////////////////////////////////////////////////////////////////////////////


function alpha = kfda(x, c, ker)
  
  class_label = unique(c);
  
  x1 = x(:, (c==class_label(1)));
  x2 = x(:, (c==class_label(2)));
  
  num = size(x, 2);
  num1 = size(x1, 2);
  num2 = size(x2, 2);
  
  M1 = sum(createkernel(x, x1, ker), 'c')/num1;
  M2 = sum(createkernel(x, x2, ker), 'c')/num2;
  M = (M1-M2)*(M1-M2)';
  
  K1 = createkernel(x, x1, ker);
  K2 = createkernel(x, x2, ker);
  l1 = ones(num1, num1)/num1;
  l2 = ones(num2, num2)/num2;
  N = K1*(eye(num1, num1)-l1)*K1'+K2*(eye(num2, num2)-l2)*K2';
  u = 100;               ////////////////////////////////////
  Nu = N + u*eye(num, num);
  
  [eigvec, eigval] = spec(inv(Nu)*M);
  [eigval, idx] = sort(diag(eigval));
  eigvec = eigvec(:, idx);
  alpha = eigvec(:, 1);
  
endfunction

