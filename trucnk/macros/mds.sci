///////////////////////////////////////////////////////////////////////////////
// Author:   Jia Wu
// Date:     Feb. 2010
// Description:  multidimensional scaling (MDS)
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
//       train_samples  -  data matrix of size dim*num; each column is a sample
//       dimension      -  output dimension with default value 2
//       criterion      -  criterion function with default 'ee'; the criterion can be:
//                         'ee'   emphasize the largest errors
//                         'ff'   emphasize the largest fractional errors
//                         'ef'   emphasize the largest product of error and fractional error
//       rate           -  convergence rate with default value 0.1
//       eps            -  accuracy with default value 0.01
//       max_iter       -  maximal number of iteration with default value 1000
//
//  Output:
//        new_samples   -  new samples after multidimensional scaling
///////////////////////////////////////////////////////////////////////////////

function new_samples = mds(train_samples, dimension, criterion, rate, eps, max_iter)
  
  //check arguments
  if argn(2)==1,
    dimension = 2;
    criterion = 'ee';
    rate = 0.1;
    eps = 0.01;
    max_iter = 1000;
  elseif argn(2)<3,
    criterion = 'ee';
    rate = 0.1;
    eps = 0.01;
    max_iter = 1000;  
  elseif argn(2)<4,
    rate = 0.1;
    eps = 0.01;
    max_iter = 1000;
  elseif argn(2)<5,
    eps = 0.01;
    max_iter = 1000;
  elseif argn(2)<6,
    max_iter = 1000;
  end
    
  [dim, num] = size(train_samples);
  
  if dimension>dim,
    error('the number of output dimension should be smaller than that of the original dimension');
  end
  
  x = train_samples;
  //initialize 'y' using principal component analysis
  y = pca2(x, dimension);
  
  //compute distance in original space
  ex = ones(1,1,num).*.x;
  delta_tmp = squeeze(sum((ex-permute(ex,[1 3 2])).^2,1));
  delta = delta_tmp(:,:);
  
  grad = zeros(y);
  delta_grad = %inf;
  iter = 0;
  
  while delta_grad>eps,
    
    //compute distances in the new space
    ey = ones(1,1,num).*.y;
    d_tmp = squeeze(sum((ey-permute(ey,[1 3 2])).^2,1));    
    d = d_tmp(:,:);
    
    for k=1:num,      
      //compute gradients
      select criterion,
        case 'ee'
          idx = [1:k-1, k+1:num];
          yk = y(:,k)*ones(1,num-1);
          yj = y(:,idx);
          d_kj = d(k,idx);
          delta_kj = delta(k,idx);
          grad(:,k) = (2/sum(triu(delta,1).^2))*sum((ones(dimension,1)*((d_kj-delta_kj)./d_kj)).*(yk-yj), 'c');
        case 'ff'
          idx = [1:k-1, k+1:num];
          yk = y(:,k)*ones(1,num-1);
          yj = y(:,idx);
          d_kj = d(k,idx);
          delta_kj = delta(k,idx);
          grad(:,k) = 2*sum((ones(dimension,1)*((d_kj-delta_kj)./(delta_kj.^2)./d_kj)).*(yk-yj), 'c');
        case 'ef'
          idx = [1:k-1, k+1:num];
          yk = y(:,k)*ones(1,num-1);
          yj = y(:,idx);
          d_kj = d(k,idx);
          delta_kj = delta(k,idx);
          grad(:,k) = (2/sum(triu(delta,1)))*sum((ones(dimension,1)*((d_kj-delta_kj)./delta_kj./d_kj)).*(yk-yj), 'c');
        else
          error('the criterion is unknown or not supported here.');
      end  
    end  
    
    old = delta_grad;
    delta_grad = sum(abs(grad));
    
    if delta_grad>3*old,
      break;
    end
        
    //update 'y'
    y = y-rate*grad;
    
    iter = iter+1;
    if iter>max_iter,
      break;
    end
    
  end

  new_samples = y;
  
endfunction

