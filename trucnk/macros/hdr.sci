///////////////////////////////////////////////////////////////////////////////
// Author:   Jia Wu
// Date:     Feb. 2010
// Description:  hierarchical dimensionality reduction (hdr)
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
//       train_samples  -  data matrix of size dim*num; each column is a sample
//       dimension      -  required output dimension 
//
//  Output:
//       new_samples    -  new data matrix of size dimension*num after dimensionality reduction
///////////////////////////////////////////////////////////////////////////////

function new_samples = hdr(train_samples, dimension)
  
  [dim, num] = size(train_samples);
  
  if dimension>=dim,
    error('the output dimension should be smaller than the original dimension of the samples');
  end
  
  d = dim;
  nd = dimension;
  
  while d~=nd,
    
    //compute correlation matrix
    sigma = train_samples*train_samples';
    mat1 = sqrt(diag(sigma)*ones(1,d));
    mat2 = sqrt(ones(d,1)*diag(sigma)');
    R_temp = sigma./(mat1.*mat2);    
    Ru = triu(R_temp, 1);
    
    //find most correlated features(dimensions)
    [val, idx] = max(Ru);
    
    //merge the two features(dimensions) by averaging them
    train_samples(idx(1), :) = (train_samples(idx(1), :)+train_samples(idx(2), :))/2
    //delete dimenson idx(2)
    train_samples = train_samples([1:idx(2)-1, idx(2)+1:d], :);
    
    d = size(train_samples, 1);
    
  end
  
  new_samples = train_samples;  
  
endfunction

