///////////////////////////////////////////////////////////////////////////////
// Author:   Jia Wu
// Date:     July 2010
// Description:  classification using a Parzen probabilistic neural network 
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
//          net             - a trained Parzen probabilistic neural network
//          test_samples    - dim*numt data matrix; each column is a data point
//          sigma           - window width of the transfer function(default 2) 
//
//  Output:
//          test_labels     - 1*numt vector of labels of each test sample
///////////////////////////////////////////////////////////////////////////////

function test_labels = pnn_classify(net, test_samples, sigma)
  
  if argn(2)<3,
    sigma = 2;
  end
  
  //augment 1 dimension to test samples
  [dimte, numte] = size(test_samples);
  test_samples(dimte+1, :) = ones(1, numte);
  dimte = dimte+1;
  
  dimtr = size(net.weight, 1);
  if dimtr~=dimte,
    error("dimension of test samples should be equal to the dimension of training samples.");
  end
  
  //normalize ???
  //tmp = sqrt(sum(test_samples.^2, 'r'));
  //test_samples = test_samples./(ones(dimte,1)*tmp);
  
  //compute the activation values
  z = net.weight'*test_samples;
  
  //transfer function
  p = transfer_func(z, sigma);
  
  //compute the discriminant 
  numc = length(net.class);
  g = zeros(numc, numte);
  for i=1:numc,
    g(i, :) = sum(p(net.classidx(i).entries, :), 'r');
  end
  
  //classify test samples according to the maximum discriminant value
  [val, ind] = max(g, 'r');
  test_labels = net.class(ind);  
  
endfunction

////////////transfer function////////////
function pr = transfer_func(z, sigma)
  
  pr = exp((z-1)./sigma.^2);
  
endfunction