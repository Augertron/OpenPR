///////////////////////////////////////////////////////////////////////////////
// Author:   Jia WU
// Date:     Nov. 2010
// Description:  Fixed-Increment Single-Sample Perceptron Criterion Function
//               (for two-category cases) 
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
//        train_samples - training data matrix of size dim*num_tr
//                        each column is a data point
//        train_labels  - 1*num_tr vector of labels for the training samples 
//        param         - parameters for the criterion: the maximum number
//                        of iterations with default value 4000
//        test_samples  - test data matrix of size dim*num_te;
//                        each column is a data point
//
//  Output:
//        a             - perceptron weights (weight vector)
//        test_labels   - predicted labels for the test samples
///////////////////////////////////////////////////////////////////////////////

function [a, test_labels] = perceptron_fis(train_samples, train_labels, param, test_samples)

    if ((argn(2)<4) & (argn(1)==2)),
        error('No input of test samples.');
    end
    
    if argn(2)<3,
    	iterm = 4000;
    else
        iterm = param(1);
    end       
    
    [dim, num_tr] = size(train_samples);

    if (length(unique(train_labels))>2),
        error('This function is for two-category cases.');
    end
    
    if (length(train_labels)~=num_tr),
        error('Number of training samples and training labels must be the same.');
    end
    
    iter = 0;
    
    //augment
    dim = dim+1;
    train_samples(dim,:) = ones(1,num_tr); 
    y = train_samples;   
    
    //"normalization"
    class = unique(train_labels);
    idx1 = find(train_labels==class(1));
    idx2 = find(train_labels==class(2));
    y(:,idx2) = -y(:,idx2);    
    
    //initialize the weights
    a = sum(y, 'c');
    
    Y_k = [1];
    
    while ((iter<iterm)&(~isempty(Y_k))),
    	
    	iter = iter+1;
    	
    	idx = ceil(rand(1)*num_tr);  //randomly choose one sample    	
    	if (a'*y(:,idx)<0),          //misclassified 
    		a = a+y(:,idx);          //updated
    	end
    	
    	Y_k = find(a'*y<0);
    	
    end    
        
    //classify
    if argn(2)==4,
    
        test_samples(dim,:) = ones(1,size(test_samples,2));  //augment
    
        if (size(test_samples,1)~=dim),
            error('Test samples and training samples should have the dimension.')
        end   
          
        test_labels = (a'*test_samples)>0;
        idx1 = find(test_labels==%T);
        idx2 = find(test_labels==%F);
        test_labels(:,idx1) = class(1);
        test_labels(:,idx2) = class(2);
        
    end
    
endfunction
