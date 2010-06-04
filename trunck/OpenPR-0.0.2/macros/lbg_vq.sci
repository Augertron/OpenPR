///////////////////////////////////////////////////////////////////////////////
// Author:   Jia Wu
// Date:     April 2010
// Description:  Linde-Buzo-Gray Vector Quantization Algorithm (LBG Design Algorithm)
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
//          training_vec   - k*M data matrix; each column is a data point
//          codevec_num    - number of code vectors
//
//  Output:
//          code_vec       - k*code_vec matrix; each column is a code vector
//          labels         - 1*M vector assigning each data point to a code vector
//          Q              - 
///////////////////////////////////////////////////////////////////////////////

function [code_vec,labels,Q] = lbg_vq(training_vec,codevec_num)

eps = 0.01;
//maximal number of iterations for each codebook size
//Max_iter = 30;

[k,M] = size(training_vec);
//number of codevectors
N = 1;

//initial codevector
code_vec = mean(training_vec,'c');

labels = zeros(1,M);

//initial Q
Q = code_vec*ones(1,M);

//average squared-error distortion
Dave = (1/(M*k))*sum((training_vec-Q).^2);

while 1

  if N<codevec_num then
    //splitting  
    code_vec = [code_vec*(1+eps),code_vec*(1-eps)];
    N = 2*N;
  elseif N>codevec_num then
    //delete one codevcector 
    //delete the codevector which represents the minimal number of training vectors
    num = zeros(1,N);
    for i = 1:N
      num(i) = length(find(labels==i));
    end;
    [val del] = min(num);
    if del(1)==1 then
      code_vec = code_vec(:,[2:$]);
    elseif del(1)==N then
      code_vec = code_vec(:,[1:$-1]);
    else
      code_vec = code_vec(:,[1:del(1)-1,del(1)+1:N]);    
    end            
    N = N-1;
  elseif N==codevec_num then
    break;
  end;

  //iteration, finding the best codevectors
  //number of iterations for each codebook size
  iter = 0;
  //distortion change between each iteration 
  delta = 1000;
  while delta>eps
    iter = iter+1;
    old_Dave = Dave;
  
    //calculate distance between training_vec and code_vec
    //and find the code_vec each training_vec belongs to
    temp = ones(1,1,N).*.training_vec;
    mat1 = permute(temp,[1,3,2]);
    mat2 = ones(1,1,M).*.code_vec;
    //N*M matrix
    dist_tmp = squeeze(sum((mat1-mat2).^2,1));
    dist = dist_tmp(:,:);
  
    //calculate labels and Q
    [val,idx] = min(dist,'r');
    labels = idx;
  
    //update codevectors
    for j = 1:N
      ind = find(labels==j);
    
      if length(ind)~=0 then
        code_vec(:,j) = mean((training_vec(:,ind)),'c');
      else //use a random vector to substitue the codevector
       r = randperm(M);
       code_vec(:,j) = training_vec(:,r(1));
      end;
    end;
  
    Q = code_vec(:,idx);

    //calculate the average error
    Dave = (1/(M*k))*sum((training_vec-Q).^2);
  
    delta = (old_Dave-Dave)/old_Dave;  

  end;

end;

endfunction
