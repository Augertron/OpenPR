///////////////////////////////////////////////////////////////////////////////
// Author:   Jia Wu
// Date:     Jan. 2010
// Description:  build a classification and regression tree
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
//       train_samples - data matrix of size dim*num; each column is a data point
//       train_labels  - class labels of training samples
//       impurity_type - impurity type for splitting; it can be Entropy, Gini, 
//                       or Misclassification
//
//  Output:
//       cart          - the trained classification and regression tree; 
//                       it is a struct variable
///////////////////////////////////////////////////////////////////////////////

function cart = buildcart(train_samples, train_labels, impurity_type) //, stop_prune_crit)

	if (length(unique(train_labels))==1),
		cart.split = train_labels(1);
		cart.right = [];
		cart.left = [];
		return;
	end
	
	ut = unique(train_labels);
	for i = 1:length(ut),
		ht(i) = length(find(labels == ut(i)));
	end
	
	if ((size(train_samples,2)==1) | (sum(ht<size(train_samples,2)))==length(ut)-1)   //stop splitting   (using Chi-square test for early stopping?)
		cart.right = [];
		cart.left = [];	
		
		[l_val, l_loc] = max(ht);
		label = ut(l_loc)
		cart.split = label;   	
	else    //split
		[dim, num] = size(train_samples);
		iN = ones(1, dim);
		split_feature = zeros(1, dim);
		threshold = zeros(1, dim);
		
		////////////////////   how to use function "optim"...
		for i = 1:dim,
			for j = 1:num,
				tmp = split(train_samples(i,j), train_samples, train_labels, i, impurity_type);
				if tmp < iN(i),
					iN(i) = tmp;
					threshold(i) = train_samples(i, j);
				end
			end
		end
		
		[val, f_dim] = min(iN);
		
		//split the node
		cart.split = [f_dim, threshold(f_dim)];
		
		right = find(train_samples(f_dim, :) > threshold(f_dim));
		left = find(train_samples(f_dim, :) <= threshold(f_dim));
		
		if isempty(right) | isempty(left),
			cart.right = [];
			cart.left = [];
			
			[l_val, l_loc] = max(ht);
			label = ut(l_loc)
			cart.split = label;   //record the label of the leaf
		else
			//continue splitting
			cart.right = buildcart(train_samples(:, right), train_labels(:, right), impurity_type) //, stop_prune_crit);
			cart.left = buildcart(train_samples(:, left), train_labels(:, left), impurity_type) //, stop_prune_crit)
		end
	end

endfunction


//split criterion
function iN = split(threshold, samples, labels, idx_f, stype)
	
	class = unique(labels);
	
	for i = 1:length(class),
		sub = find(labels == class(i));
		Pr(i) = length(find(samples(idx_f, sub) > threshold))/length(sub);
		Pl(i) = length(find(samples(idx_f, sub) <= threshold))/length(sub);
	end
	
	select stype,
		case 'Entropy'
			iN_r = sum(-Pr.*log(Pr+2.2204e-016)/log(2));
			iN_l = sum(-Pl.*log(Pl+2.2204e-016)/log(2));
		case 'Gini'
			iN_r = 1 - sum(Pr.^2);
			iN_l = 1 - sum(Pl.^2);
		case 'Missclassification'
			iN_r = 1 - max(Pr);
			iN_l = 1 - max(Pl);		
		otherwise
			error('wrong split type');
	end
	
	ra = length(find(samples(idx_f, :) > threshold))/length(labels);
	iN = -ra*iN_r-(1-ra)*iN_l;

endfunction


//stop criterion


//prune criterion
