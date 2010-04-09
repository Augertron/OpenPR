///////////////////////////////////////////////////////////////////////////////
// Author:  Jia Wu
// Version: 0.1
// Date: Nov 2009
// Description: Multiple Linear Discriminant Analysis(LDA)
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
//        x       - Data matrix of size dim*num. Each column is a data point.
//        c 	  - Class label vector of size 1*num or num*1.
//		
//  Output:
//        w       - Weight matrix of size dim*(class_num-1).
///////////////////////////////////////////////////////////////////////////////

function w = mlda(x, c)
	
	[dim, num] = size(x);
	
	class_label = unique(c);
	class_num = length(class_label);
	
	if(dim < class_num-1),
	    error('Number of class must be no bigger than the data dimension.');
    end
  
	if dim > num-class_num,	//in case of Sw is sigular
		[new_x, mn, eig_val, eig_vec] = pca2(x, num-class_num);
		w0 = lda(new_x, c);
		
		w = eig_vec*w0;
	else	
		m0 = mean(x, 2);
		Sw = zeros(dim, dim);
		Sb = zeros(dim, dim);
		
		for i=1:class_num,
			mat = x(:, (c==class_label(i)));
			n = size(mat, 2);
			m = mean(mat, 'c');
			S = (mat-m*ones(1,n))*(mat-m*ones(1,n))';
			
			Sw = Sw + S;
			Sb = Sb + n*(m-m0)*(m-m0)';
		end
		
        [al, be, w] = spec(Sb, Sw);
        lambda = al./be;
		[lambda, idx] = sort(diag(lambda));
		w = w(:, idx);
		
		w = w(:, 1:(class_num-1));
        end

endfunction
