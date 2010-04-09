///////////////////////////////////////////////////////////////////////////////
// Author:   Jia Wu
// Date:     Feb. 2010
// Description:  competitve learning clustering
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
//        train_samples - data matrix of size dim*num; each column is a data point
//        c             - number of clusters
//        eta           - learning rate with default value 0.01
//        alpha         - decay coefficient of learning rate with default value 0.99
//        max_iter      - maximal number of iterations with default value 1000
//        eps           - threshold for change in weight vector with default value 1e-5
//
//  Output:
//        centers       - cluster centers
//        labels        - cluster indices for each training sample point
//        W             - weight vectors
///////////////////////////////////////////////////////////////////////////////

function [centers, labels, W] = competitive_learning(train_samples, c, eta, alpha, max_iter, eps)
  
  if argn(2)<3,
    eta = 0.01;
    alpha = 0.99;
    max_iter = 1000;
    eps = 1e-5;
  elseif argn(2)<4,
    alpha = 0.99;
    max_iter = 1000;
    eps = 1e-5;
  elseif argn(2)<5,
    max_iter = 1000;
    eps = 1e-5;
  elseif argn(2)<6,
    eps = 1e-5;
  end
  
  [dim, num] = size(train_samples);
  
  //augment and normalize samples
  x = [ones(1, num); train_samples];
  x = x./(ones(dim+1, 1)*sqrt(sum(x.^2, 'r')));
  
  //initialize weight vectors
  [val, idx] = sort(rand(1, num));
  sel = idx(1:c);
  W = x(:, sel);  
  
  //learning
  iter = 0;
  change = 1000;
  while change > eps,
    
    iter = iter+1;
    change = 0;
    
    //randomly reorder the samples
    [val, idx] = sort(rand(1, num));
    
    for i = 1:num
      J = W'*x(:,idx(i));
      [val, j] = max(J);
      
      old_w = W(:,j);      
      //weight update
      W(:,j) = W(:,j)+eta*x(:,idx(i));
      //weight normalization
      W(:,j) = W(:,j)/sqrt(sum(W(:,j).^2, 'r'));
      
      change = change+sum(abs(W(:,j)-old_w));      
    end
        
    if iter > max_iter,
      break;
    end
  
    change = change/num;
    eta = eta*alpha;
    
  end
    
  //compute cluster center and assign sample labels
  centers = W(2:$, :);  
  res = W'*x;
  [val, labels] = max(res, 'r');
 
endfunction

