///////////////////////////////////////////////////////////////////////////////
// Author:   Jia Wu
// Date:     Feb. 2010
// Description:  basic iterative minimum-squared-error clustering
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

function [centers, labels] = mseclustering(train_samples, cluster_num) //, attempt_num)
  
//  if argn(2)<3,
//    attempt_num = 1;
//  end
  
  num = size(train_samples, 2);
  attempt_num = num;    //number of attempts
  
  //initialize cluster centers
  [val, idx] = sort(rand(1, num));
  centers = train_samples(:, idx(1:cluster_num));
  
  //initial partition of the samples according to nearest distance
  dist = zeros(cluster_num, num);
  for i=1:cluster_num,
    dist(i, :) = sum(((centers(:, i)*ones(1, num))-train_samples).^2, 'r');
  end
  [val labels] = min(dist, 'r');
  
  n = zeros(1, cluster_num);
  for i=1:cluster_num,
    n(i) = length(find(labels==i));
  end

  //sum-of-squared-error
  Je = 0;
  
  //iterative minimum-squared-error clustering
  while attempt_num,
    
    old_Je = Je;
        
    //randomly select a sample
    [val, idx] = sort(rand(1, num));
    x = train_samples(:, idx(1));    
    
    //classify sample x
    dist = sum((x*ones(1, cluster_num)-centers).^2, 'r');
    [val, idx1] = find(dist==min(dist));
    lables(idx(1)) = idx1;
    
    ro = zeros(1, cluster_num);
    
    if n(idx1)~=1,
      ro = (n./(n+1)).*dist;
      ro(idx1) = ro(idx1)*((n(idx1)+1)/(n(idx1)-1));
 
      [val, idx2] = find(ro==min(ro));
      
      //transfer samples x and recompute Je and cluster centers
      if idx2~=idx1,
        labels(idx(1)) = idx2;
        n(idx2) = n(idx2)+1;
        n(idx1) = n(idx1)-1;
        
        for i=1:cluster_num,
          centers(:, i) = mean(train_samples(:,find(labels==i)), 'c');
          
		  J(i) = sum((train_samples(:,find(labels==i))-centers(:,i)*ones(1,length(find(labels==i)))).^2);
        end      
        Je = sum(J);
      end
    end  
  
//  if Je==old_Je,
//    attempt_num = attempt_num-1;
//  end
  
    if Je~=old_Je,
      attempt_num = attermpt_num-1;
    else
      break;
    end
    
  end
  
endfunction

