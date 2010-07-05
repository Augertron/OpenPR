///////////////////////////////////////////////////////////////////////////////
// Author:   Jia Wu
// Date:     June 2010
// Description:  Mean-Shift Clustering
// Reference:    Fukunaga, K.; Hostetler, L.; , "The estimation of the gradient
//               of a density function, with applications in pattern recognition,"
//               Information Theory, IEEE Transactions on , vol.21, no.1, 
//               pp. 32- 40, Jan 1975
//
// Copyright (C) 2009-2010 OpenPR
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
//          samples       - dim*num data matrix; each column is a data point
//          h             - parameter for bandwidth
//
//  Output:
//          cluster_centers       - cluster centers
//          sample_id     - label for each sample to which cluster it belongs to
///////////////////////////////////////////////////////////////////////////////

function [cluster_centers, sample_id] = msclustering(samples, h)
  
  [dim, num] = size(samples);
  
  //threshold for convergence
  esp = 1e-3*h;
 
  used_flag = zeros(1, num);
  clusters_count = zeros(1, num);
  cluster_num = 0;
  cluster_centers = [];
  
  num_unused = num;
  sp_unused = 1:num;
  
  while num_unused
    
    //randomly pick a sample as an initialization
    idx = randperm(num_unused);
    X = samples(:, sp_unused(idx(1)));

    onecluster_count = zeros(1, num);
    
    while 1, //iterative loop finding the mode
      
      //find samples inside the window
      dist = sum((samples-X*ones(1,num)).^2, 'r');
      sp_idx = find(dist<=h^2);
      onecluster_count(sp_idx) = onecluster_count(sp_idx)+1;
      used_flag(sp_idx) = 1;      
      
      old_X = X;
      //formulae: 
      //Xi+1 = 1/k*ΣX,  X∈Sh(Xi);Sh(X)={Y:(Y-X)'(Y-X)<=h^2}
      X = mean(samples(:,sp_idx), 'c');
      
      //if convergence
      if sum((X-old_X).^2)<esp,
        
        //calculate distances between X and previous cluster centers
        c_dist = sum((cluster_centers-X*ones(1,cluster_num)).^2, 'r');
        merge_idx = find(c_dist<=h/2);
        
        if length(merge_idx)==0,
          //add a new cluster center
          cluster_num = cluster_num+1;
          cluster_centers(:,cluster_num) = X; 
          cluster_count(cluster_num, :) = onecluster_count;
        else
          //merge the closest two
          [temp_val temp_idx] = min(c_dist(merge_idx));
          cluster_centers(:,merge_idx(temp_idx)) = (X+cluster_centers(:,merge_idx(temp_idx)))/2;
          cluster_count(merge_idx(temp_idx), :) = cluster_count(merge_idx(temp_idx), :)+onecluster_count;
        end        
        
        break;
        
      end
                  
    end    
    
    //find unused samples
    sp_unused = find(used_flag==0);
    num_unused = length(sp_unused);
    
  end
  
  [val, sample_id] = max(cluster_count, 'r');
 
 
endfunction