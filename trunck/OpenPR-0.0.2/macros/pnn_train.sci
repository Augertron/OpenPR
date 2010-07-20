///////////////////////////////////////////////////////////////////////////////
// Author:   Jia Wu
// Date:     July 2010
// Description:  create a Parzen probabilistic neural network
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
//
//  Output:
//          net             - Parzen probabilistic neural network 
///////////////////////////////////////////////////////////////////////////////

function net = pnn_train(train_samples, train_labels) 
  
  [dim, nums] = size(train_samples);
  numl = length(train_labels);
  
  if nums~=numl,
    error("the number of samples and the number of labels must be equal.");
  end  
  
  //augment 1 dimension to the training samples
  train_samples(dim+1, :) = ones(1, nums);
  dim = dim+1;
  
  class = unique(train_labels);
  numc = length(class);
  
  //normalize
  x_tmp = sqrt(sum(train_samples.^2, 'r'));
  x = train_samples./(ones(dim,1)*x_tmp);
  
  //create the network
  net.weight = x;
  net.class = class;
  net.classidx = cell(1, numc);
  for i=1:numc,
    net.classidx(i).entries = find(train_labels==class(i));
  end
    
endfunction