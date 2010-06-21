///////////////////////////////////////////////////////////////////////////////
// Author:   Jia Wu
// Date:     June 2010
// Description:  use least-squares algorithm to classify 
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
//         train_samples      - dim*num data matrix; 
//                              each column is a data sapmle, each row is a feature
//         train_labels       - 1*num vector of labels for each data sample
//         test_samples       - dim*test_num data matrix
//
//  Output:
//         test_labels        - 1*test_num vector of predicted labels for test_samples
//         a                  - weighted vector
//
///////////////////////////////////////////////////////////////////////////////

function [test_labels, a] = least_squares(train_samples, train_labels, test_samples)
  
  [dim num] = size(train_samples);
  [test_dim test_num] = size(test_samples);
  
  if test_dim~=dim,
    error('the dimension of test samples must equal to that of training samples');
  end
  
  //augment
  train_samples(dim+1, :) = ones(1, num);
  test_samples(dim+1, :) = ones(1, test_num); 
  
  //calculate the weighted vector: a
  //  Ya = b ---> Y'Ya = Y'b ---> a = (Y'Y)^(-1)Y'b
  //  (Y: each row is data sample)
  a = (inv(train_samples*train_samples')*train_samples*train_labels')';
  
  //classify test_samples
  test_labels = a*test_samples;
  test_labels = floor(test_labels);
  
endfunction
