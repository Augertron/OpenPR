///////////////////////////////////////////////////////////////////////////////
// Author:  Jia Wu
// Version: 0.1
// Date: Nov 2009
// Description: Independent Component Analysis(ICA)
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
//        mixture       - observed signals; each column is a data point
//        K             - number of source components
//		
//  Output:
//		  source        - source signals
//        W             - weight matrix
//        Aw            - whitening matrix
///////////////////////////////////////////////////////////////////////////////


function [source, W, Aw] = ica(mixture, K)

	[r, c] = size(mixture);

	//scale 
	mixture_org = mixture;
	scale = max(abs(mixture));
	mixture = mixture./scale;

	//use SVD to reduce dimension
	if(argn(2) < 2),
		K = r;
	end 
	 
	if(K > r),
		error('number of source components cannot be larger than the number of observed signals')
	end

	if((K > 0) & (K < r)),
		if(r < c),
			[V, S, U] = svd(mixture', 'e');
		else
			[U, S, V] = svd(mixture, 'e');
		end
		
		mixture = S(1:K, 1:K)*V(:, 1:K)';
	end

	r1 = size(mixture, 1); 

	//Whitening transformation
	[mixture, Aw] = whitening(mixture);

	//calculate weight matrix
	W = rand(r1, r1);
	iter = 1;
	eta = 0.01;
	while(iter < 1000),
		iter = iter + 1;
		y = W*mixture;
		dW = (eye(r1,r1)-1/c*(y.^3)*y')*W;                //

		if(max(dW) > 1000),
			disp('Algorithm diverged.');
			break;
		end

		W = W + eta*dW;                                   //

		update = max(abs(dW)); 
		if(update < eta),
			//disp('Algorithm converged.');
			break;
		end
	end
		
	[val, idx] = sort(sum(abs(W), 'c'));
	W = W(idx, :);

	//calculate the source signals and the weight matrix
	source = W*mixture;
	W = W*Aw';           

endfunction
