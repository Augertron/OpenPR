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
//        stop         - a scalar in (0, 1). Stop iterations if improved rate is
//                       less than this value.
//  Output:
//        labels       - labels of the input data.
//        centroids    - cluster centroids
///////////////////////////////////////////////////////////////////////////////

function [labels, centroids] = kmeans_sci(data, k, stop)

	if k<1,
		error('k should be a positive integer.');
	end
	
	[dim, num] = size(data);
	labels = zeros(1, num);
	dist = zeros(k, num);
	
	//initialize means
	[val, idx] = sort(rand(1, num));
	sel = idx(1:k);
	m = data(:, sel);

	if k == 1,
		m = mean(data, 'c');
		labels = ones(1, num);
		return;
	end
	
	rate = 1e6;
	ratio = 1e6;
	while(ratio > stop),
	
		for i=1:k,
			dist(i,:) = sum((data-m(:,i)*ones(1,num)).^2, 'r');
		end
		
		[minimun, labels] = min(dist, 'r');
		
        //recompute means
		for i=1:k,
		    m(:,i) = mean(data(:, find(labels==i)), 'c');
        end
      
        old_rate = rate;
      		rate = mean(minimun);
      		ratio = 1-rate/old_rate;
      		
     end

     centroids = m;
     
endfunction
