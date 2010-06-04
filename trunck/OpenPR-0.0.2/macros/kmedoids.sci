///////////////////////////////////////////////////////////////////////////////
// Author:   Jia Wu
// Date:     May 2010
// Description:  basic k-medoids clustering algorithm
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
//          samples - dim*snum data matrix; each column is a data point
//          knum    - number of the medoids
//
//  Output:
//          medoids - dim*knum matrix; each column is a cluster center
//          labels  - 1*snum vector assigning each data point to a cluster
///////////////////////////////////////////////////////////////////////////////

function [medoids, labels] = kmedoids(samples, knum)
   
[dim, snum] = size(samples);

medoids = zeros(dim, knum);

labels = zeros(1, snum);      

//initialize medoids; randomly select 'knum' of sample points as cluster centers
idx = randperm(snum);
medoids = samples(:, idx(1:knum));
new_medoids = samples(:, idx(1:knum)+1);

while max(abs(new_medoids-medoids))>0.01,
  
  medoids = new_medoids;
  
  //calculate distance using the Euclidean distance;it could be other distance metric like Manhattan distance and Minkowski distance
  dist_matrix = zeros(knum, snum);        
  for i=1:knum,
    dist_matrix(i,:) = sum((samples-medoids(:,i)*ones(1, snum)).^2, 1);
  end
  
  [temp, idx] = min(dist_matrix, 'r');
  labels = idx;
  
  //swap cluster centers and sample points; select the configuration with the lowest cost and change the medoids
  for i=1:knum,
    tmp = find(idx==i);
    if ~isempty(tmp),
      cluster_samples = samples(:, tmp);
      num = size(cluster_samples, 2);
      samples_matrices = ones(1,1,num).*.cluster_samples;
      centers_matrices = permute(samples_matrices, [1 3 2]);
      cost = sum(sum((samples_matrices-centers_matrices).^2,1),2);
      [val, cen_idx] = min(cost, 3);
      new_medoids(:, i) = cluster_samples(:, cen_idx(:,:,:));
    else
      new_medoids(:,i) = samples(:, floor(rand(1)*snum)+1);
    end  
  end
  
end
            
endfunction