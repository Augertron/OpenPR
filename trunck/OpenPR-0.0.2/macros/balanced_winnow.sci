///////////////////////////////////////////////////////////////////////////////
// Author:   Jia WU
// Date:     Nov. 2010
// Description:  Balanced Winnow Algorithm (for two-category cases) 
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
//        param         - parameters for the criterion, including the maximum
//                        number of iterations, the convergence rate and the
//                        promotion parameter: [iterm,eta,alpha] 
//                        The default value is [1000,0.01,2]
//        test_samples  - test data matrix of size dim*num_te
//                        each column is a data point
//
//  Output:
//        a_plus        - positive weight vector
//        a_minus       - negative weight vector
//        test_labels   - predicted labels for the test samples
///////////////////////////////////////////////////////////////////////////////

function [a_plus, a_minus, test_labels] = balanced_winnow(train_samples, train_labels, param, test_samples)

    if ((argn(2)<4) & (argn(1)==3)),
        error('No input of test samples.');
    end
    
    if argn(2)<3,
    	iterm = 1000;
    	eta = 0.01;
    	alpha = 2;
    else
        iterm = param(1);
        eta = param(2);
        alpha = param(3);
    end    
    
    [dim, num_tr] = size(train_samples);

    if (length(unique(train_labels))>2),
        error('This function is for two-category cases.');
    end
    
    if (length(train_labels)~=num_tr),
        error('Number of training samples and training labels must be the same.');
    end
    
    //augment
    dim = dim+1;
    train_samples(dim,:) = ones(1,num_tr); 
    y = train_samples;   
    z = train_labels;
    
    //"normalization"
    class = unique(train_labels);
    idx1 = find(train_labels==class(1));
    idx2 = find(train_labels==class(2));
    //y(:,idx2) = -y(:,idx2);    
    
    z(:,idx1) = 1;
    z(:,idx2) = -1;
    
    //initialize the weights
    //a_plus = 1/(2*dim)*ones(dim,1);
    //a_minus = -a_plus;
    a_plus = sum(y,'c');
    a_minus = -sum(y,'c');

    iter = 0;
        
    while (iter<iterm),
    	
        iter = iter+1;
        
        for i=1:num_tr,
          
          if (sign(a_plus'*y(:,i)-a_minus'*y(:,i))~=z(i)),  //misclassified sample 
            //update weights
            if z(i)==1,
              a_plus = alpha.^(eta*y(:,i)).*a_plus;
              a_minus = alpha.^-(eta*y(:,i)).*a_minus;
            else
              a_plus = alpha.^-(eta*y(:,i)).*a_plus;
              a_minus = alpha.^(eta*y(:,i)).*a_minus;
            end
            
            Y_k = find((a_plus'*y-a_minus'*y)~=z);
            if isempty(Y_k),
              break;
            end
            
          end
          
        end        
    
    end    
        
    //classify
    if argn(2)==4,
    
        if (size(test_samples,1)~=dim-1),
            error('Test samples and training samples should have the dimension.')
        end   
        
        test_samples(dim,:) = ones(1,size(test_samples,2));  //augment
        
        a = (a_plus+a_minus)/2; 
          
        test_labels = (a'*test_samples)<0.5;
        idx1 = find(test_labels==%T);
        idx2 = find(test_labels==%F);
        test_labels(:,idx1) = class(1);
        test_labels(:,idx2) = class(2)
        
    end
    
endfunction
