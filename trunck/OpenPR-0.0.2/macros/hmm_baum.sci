///////////////////////////////////////////////////////////////////////////////
// Author:   Jia Wu
// Date:     August 2010
// Description:  Baum-Welch algorithm for estimating the parameters of a
//               hidden markov model given an observation sequence
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
//  Iutput:
//          num_states      - number of hidden states
//          num_symbols     - number of visible states
//          ob_seq          - observation sequence in numerics
//
//  Onput:
//          tran_prob       - transition probability matrix
//          emit_prob       - emission probability matrix
//          start_prob      - initial probability vector
///////////////////////////////////////////////////////////////////////////////

function [tran_prob, emit_prob, start_prob] = hmm_baum(num_states, num_symbols, ob_seq)
  
  N = num_states;
  M = num_symbols;
  T = length(ob_seq);
  
  //initialize the transition probabilities, emission probabilites and the 
  //initialization probabilites of the model
  tran_prob = rand(N,N);
  emit_prob = rand(N,M);
  start_prob = rand(1,N);
  
  tran_prob = tran_prob./(sum(tran_prob, 'c')*ones(1,N));
  emit_prob = emit_prob./(sum(emit_prob, 'c')*ones(1,M));
  start_prob = start_prob/sum(start_prob);  
  
  [P, alpha] = hmm_forward(tran_prob, emit_prob, ob_seq, start_prob);
  [Pb, be] = hmm_backward(tran_prob, emit_prob, ob_seq, start_prob);
  ga = compute_gamma(T, N, alpha, be);
  xi = compute_xi(T, N, ob_seq, tran_prob, emit_prob, alpha, be);
      
  delta = 1;
  eps = 0.001;
  
  while delta>eps,
          
    //reestimate
    //initialization probabilities      
    start_prob = 0.001+0.999*ga(1,:);
    
    //transition probabilities
    for i=1:N,
      for j=1:N,
        tran_prob(i,j) = 0.001+0.999*(sum(xi(i,j,:))/sum(ga(1:T-1,i)));
      end
    end
    
    //emission probabilities
    for j=1:N,
      for k=1:M,
        idx = find(ob_seq==k);
        emit_prob(j,k) = 0.001+0.999*(sum(ga(idx,j))/sum(ga(:,j)));
      end
    end
    
    P_old = P;
    [P, alpha] = hmm_forward(tran_prob, emit_prob, ob_seq, start_prob);
    [Pb, be] = hmm_backward(tran_prob, emit_prob, ob_seq, start_prob);
    ga = compute_gamma(T, N, alpha, be);
    xi = compute_xi(T, N, ob_seq, tran_prob, emit_prob, alpha, be);
    
    delta = abs(P-P_old);
    
  end
    
endfunction

///////////compute gamma: gamma(t,i) = P(q(t)=Si|Ot)///////////
function ga = compute_gamma(t, n, alpha, be)
    
  ga = zeros(t, n);
  temp = alpha.*be;
  ga = temp./(sum(temp, 'c')*ones(1,n));
  
endfunction

///////////compute xi: xi(i,j,t) = P(q(t)=Si, q(t+1)=Sj|Ot)///////////
function xi = compute_xi(t, n, ob_seq, tran_prob, emit_prob, alpha, be)

  xi = zeros(n, n, t-1);
  
  for tm=1:t-1,
    tmp = 0;
    for i=1:n,
      for j=1:n,
        xi(i,j,tm) = alpha(tm,i)*tran_prob(i,j)*emit_prob(j,ob_seq(tm+1))*be(tm+1,j);
        tmp = tmp+xi(i,j,tm);
      end      
    end
    xi(:,:,tm) = xi(:,:,tm)/tmp;
  end
    
endfunction         
