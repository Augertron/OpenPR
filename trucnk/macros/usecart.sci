///////////////////////////////////////////////////////////////////////////////
// Author:   Jia Wu
// Date:     Jan. 2010
// Description:  use a trained classification and regression tree to predict 
//               the labels of given samples
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
//       samples  - input sample matrix; each column is a sample point
//       indices  - indices of selected samples
//       cart     - the trained classification and regression tree by function 'buildcart'
//
//  Output:
//       labels   - predicted labels of the input samples
///////////////////////////////////////////////////////////////////////////////

function labels = usecart(samples, indices, cart)

	num = size(samples, 2);
	
	if length(cart.split) == 1,  //leaf
		labels = zeros(1, num);
		labels(indices) = cart.split;
	else                         //node
		right_indices = indices(find(samples(cart.split(1), indices) > cart.split(2)));
		left_indices = indices(find(samples(cart.split(1), indices) <= cart.split(2)));
		
		if isempty(right_indices)
			rlabels = zeros(1, num);
		else
			rlabels = usecart(samples, right_indices, cart.right);
		end
			
		if isempty(left_indices)
			llabels = zeros(1, num);
		else
			llabels = usecart(samples, left_indices, cart.left);
		end
		
		labels = rlabels + llabels;
	end

endfunction
