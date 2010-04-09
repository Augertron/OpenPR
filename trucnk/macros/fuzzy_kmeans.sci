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

///////////////////////////////////////////////////////////////////////////////
//  Input:
//        data         - dim*num data matrix; each column is a data point.
//        k            - number of nearest neighbors.
//        b            - b > 1 is a free parameter chosen to adjust the 
//                       ¡°blending¡± of different clusters
//  Output:
//        labels       - labels of the input data.
//        centroids    - cluster centroids
///////////////////////////////////////////////////////////////////////////////

function [labels, centroids] = fuzzy_kmeans(data, k, b)
  
  if k<1,
    error('k should be a positive integer.');
  end
  
  [dim, num] = size(data);
  labels = zeros(1, num);
  dist = zeros(k, num);
  
  //initialize means
  [val, idx] = sort(rand(1, num));
  sel = idx(1:k);
  m = data(:, sel)+mean(data, 'c')*ones(1,k);
  o_m = zeros(dim, k);
  
  //initialize degrees
  p = rand(k, num);
  p = p./(ones(k,1)*sum(p,'r'));
  o_p = zeros(k, num);
  
  while((sum(abs(m-o_m))>1e-5)&(sum(abs(p-o_p))>1e-5)),
    
    o_m = m;
    o_p = p;
    
    for i=1:k,
      dist(i,:) = sum((data-m(:,i)*ones(1,num)).^2, 'r');
    end    
    
    //recompute degrees
    p = (dist.^(-1)).^(1/(b-1));
    p = p./(ones(k,1)*sum(p,'r')); 
    
    //recompute means
    p = p.^b;
    m = (data*p')./(ones(dim,1)*sum(p, 'c')');     
    
  end  
  
  [maximun, labels] = max(p, 'r');
  centroids = m;
  
endfunction
