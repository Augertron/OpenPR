//////////////////////////////////////////////////////////////////////////////
// Author:  Jia Wu
// Version: 0.1
// Date: Dec. 2009
// Description: K-Nearest Neighbors
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
//        train_data         - dim*num data matrix; each column is a data point.
//		  train_labels		 - labels of each train_data.
//        test_data          - dim*t_num data matrix; each column is a data point.
//        k                  - number of nearest neighbors.
//  Output:
//        test_labels        - labels of each test_data.
///////////////////////////////////////////////////////////////////////////////

function test_labels = knearest_sci(train_data, train_labels, test_data, k)

	num = size(train_data, 2);
	if k>num, 
		error('number of nearest neighbors must not be bigger than the number of train_data.');
	end
	
	t_num = size(test_data, 2);	            //sample number
	test_labels = zeros(1, t_num);
	
	for i=1:t_num,
		//Euclidean Distance
		dist = sum((test_data(:, i)*ones(1, num)-train_data).^2, 'r'); 
		
		[vec, idx] = sort(dist);
		labels = train_labels(idx(k+1:num)); 
		
		tmp_labels = unique(labels);
		tmp_num = length(tmp_labels);
		tmp_vec = zeros(1, tmp_num);
		for j=1:tmp_num,
			tmp_vec(j) = sum(labels==tmp_labels(j));
		end
		
		[mx, idx] = max(tmp_vec);
		test_labels(i) = tmp_labels(idx);			
	end

endfunction
