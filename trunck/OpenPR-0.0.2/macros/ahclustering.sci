///////////////////////////////////////////////////////////////////////////////
// Author:   Jia Wu
// Date:     Feb. 2010
// Description:  agglomerative hierarchical clustering
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
//        dist_type      - type of distance function with default value 'min'
//                         it could be 'min', 'max', 'avg' or 'mean'
//
//  Output:
//        centers        - centers of the formed clusters
//        labels         - labels of each trainning sample belonging to the formed clusters
///////////////////////////////////////////////////////////////////////////////

function [centers, labels] = ahclustering(train_samples, cluster_num, dist_type)
  
  if argn(2)<3,
    dist_type = 'min';
  end

  [dim, num] = size(train_samples);
  
  labels = [1:num];  
  current_num = num;    //current number of clusters
  
  //calculate distance between each pair of samples 
  mat1 = ones(1,1,num).*.train_samples;
  mat2 = permute(mat1, [1 3 2]);
  dist = sqrt(squeeze(sum((mat1-mat2).^2, 1)));
  
  while current_num>cluster_num,
    
    class = unique(labels);
    dist_ab = zeros(current_num, current_num);
    
    select dist_type,
      case 'min'
        for i=1:current_num,
          idx_a = find(labels==class(i));
          for j=1:current_num,
            idx_b = find(labels==class(j));
            
            dist_mat = dist(idx_a, idx_b);
            dist_ab(i,j) = min(dist_mat);
          end          
        end
        dist_ab = dist_ab+1e20*eye(current_num,current_num);  //exclude samples in the same cluster
        [a, b] = find(dist_ab == min(dist_ab));
        a = class(a(1));
        b = class(b(1));
      case 'max'
        for i=1:current_num,
          idx_a = find(labels==class(i));
          for j=1:current_num,
            idx_b = find(labels==class(j));
            
            dist_mat = dist(idx_a, idx_b);
            dist_ab(i,j) = max(dist_mat);
          end          
        end
        dist_ab = dist_ab.*(ones(current_num,current_num)-eye(current_num,current_num));  //exclude samples in the same cluster
        [a, b] = find(dist_ab == max(dist_ab));
        a = class(a(1));
        b = class(b(1));
      case 'avg'
        for i=1:current_num,
          idx_a = find(labels==class(i));
          for j=1:current_num,
            idx_b = find(labels==class(j));
            
            dist_mat = dist(idx_a, idx_b);
            dist_ab(i,j) = mean(dist_mat);
          end          
        end
        dist_ab = dist_ab+1e20*eye(current_num,current_num);  //exclude samples in the same cluster
        [a, b] = find(dist_ab == min(dist_ab));
        a = class(a(1));
        b = class(b(1));        
      case 'mean'   
        for i=1:current_num,
          set_a = train_samples(:, find(labels==class(i)));
          m_a = mean(set_a, 'c');
          for j = 1:current_num,
            set_b = train_samples(:, find(labels==class(j)));
            m_b = mean(set_b, 'c');
            
            dist_ab(i, j) = sqrt(sum((m_a-m_b).^2, 'r'));
          end
        end
        dist_ab = dist_ab+1e20*eye(current_num,current_num);  //exclude samples in the same cluster
        [a, b] = find(dist_ab == min(dist_ab));
        a = class(a(1));
        b = class(b(1));           
      else
        error('distance measures should be one of min, max, avg, or mean.');
    end    
    
    //merge the nearest two clusters
    labels(find(labels==b)) = a;
    current_num = length(unique(labels));         
        
  end

  //calculate cluster centers
  centers = zeros(dim, cluster_num);
  class = unique(labels);
  for i=1:cluster_num,
    idx = find(labels == class(i));
    centers(:, i) = mean(train_samples(:,idx), 'c');
  end
   
endfunction

