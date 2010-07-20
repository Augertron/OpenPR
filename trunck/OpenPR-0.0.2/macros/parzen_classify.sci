///////////////////////////////////////////////////////////////////////////////
// Author:   Jia Wu
// Date:     July 2010
// Description:  classification using simple parzen-window estimation 
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
//          train_samples   - dim*num data matrix; each column is a data point
//          train_labels    - 1*num vector of labels of each training sample
//          test_samples    - dim*numt data matrix; each column is a data point
//          h               - window width
//
//  Output:
//          test_labels     - 1*numt vector of labels of each test sample 
///////////////////////////////////////////////////////////////////////////////

function test_labels = parzen_classify(train_samples, train_labels, test_samples, h)
  
  num = size(train_samples, 2);
  numt = size(test_samples, 2);
  class = unique(train_labels);
  
  P = zeros(1, length(class));
  k = zeros(length(class), numt);
  
  for i=1:length(class),
    sample_num = length(find(train_labels==class(i)));
    P(i) = sample_num/num;    
    tmp_samples = train_samples(:, find(train_labels==class(i)));
    
    for j=1:numt,      
      dist = sqrt(sum((test_samples(:,j)*ones(1, sample_num)-tmp_samples).^2, 'r'));
      k(i,j) = sum(phi(dist./h));
    end
    
    k(i,:) = (k(i,:)/sample_num)*P(i);
  end
  
  //classification according to the maximum posterior
  [val, idx] = max(k, 'r');
  test_labels = idx;      
  
endfunction


///////window function//////////////////////////////
function k = phi(mu)
  
  k = (abs(mu)<=1/2);
  
endfunction
