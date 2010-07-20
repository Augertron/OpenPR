///////////////////////////////////////////////////////////////////////////////
// Author:   Jia Wu
// Date:     July 2010
// Description:  classification using a reduced coulomb energy network 
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
//          net             - a trained reduced coulomb energy network
//          test_samples    - dim*numt data matrix; each column is a data point
//
//  Output:
//          test_labels     - 1*numt vector of labels of each test sample
///////////////////////////////////////////////////////////////////////////////

function test_labels = rce_classify(net, test_samples)
  
  [dim, numt] = size(test_samples);
    
  if dim~=size(net.weight,1),
    error("dimension of test samples should be equal to the dimension of training samples.");
  end
  
  for i=1:numt,
    dist = sqrt(sum((test_samples(:,i)*ones(1,size(net.weight,2))-net.weight).^2, 'r'));
    idx = find(dist<net.lambda);
    
    if isempty(idx),
      test_labels(i) = 0;  //////////////
    else
      class = unique(net.label(idx));
      if length(class)==1,
        test_labels(i) = class;
      else
        test_labels(i) = 0;   //////////////
      end      
    end
    
  end   
  
  test_labels = test_labels';
endfunction