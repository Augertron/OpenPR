///////////////////////////////////////////////////////////////////////////////
// Author:   Jia Wu
// Date:     Feb. 2010
// Description:  stepwise optimal hierarchical clustering
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
//        train_samples  - data matrix of size dim*num; each column is a data point
//        cluster_num    - number of desired clusters
//
//  Output:
//        centers        - centers of the formed clusters
//        labels         - labels of each trainning sample belonging to the formed clusters
///////////////////////////////////////////////////////////////////////////////
function [centers, labels] = sohclustering(train_samples, cluster_num)
  
  [dim, num] = size(train_samples);
  
  labels = [1:num];
  current_num = num;    //current number of clusters
  
  while current_num>cluster_num,
    
    cluster = unique(labels);
    de = zeros(current_num, current_num);
    
    for i=1:current_num,
      idx_a = find(labels==cluster(i));
      m_a = mean(train_samples(:,idx_a), 'c');
      for j = 1:current_num,
        idx_b = find(labels==cluster(j));
        m_b = mean(train_samples(:,idx_b), 'c');
        
        de(i, j) = (length(idx_a)*length(idx_b)/(length(idx_a)+length(idx_b)))*sqrt(sum((m_a-m_b).^2, 'r'));
      end
    end
    
    de = de+1e20*eye(current_num, current_num);
    [a, b] = find(de==min(de));
    a = cluster(a(1));
    b = cluster(b(1));
    
    //merge clusters
    labels(find(labels==b)) = a;
    current_num = current_num-1;  //current_num = length(unique(labels))
    
  end
  
  centers = zeros(dim, cluster_num);
  cluster = unique(labels);
  for i=1:cluster_num,
    idx = find(labels == cluster(i));
    centers(:, i) = mean(train_samples(:,idx), 'c');
  end
  
endfunction

