///////////////////////////////////////////////////////////////////////////////
// Author:   Jia Wu
// Date:     July 2010
// Description:  Viterbi algorithm for computing the optimal hidden state
//               sequence given an observation sequence and the hidden markov model
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
//          tran_prob       - transition probability matrix of size N*N
//          emit_prob       - emission probability matrix of size N*M
//          ob_seq          - observation sequence in numerics
//          start_prob      - initial probability vector of size 1*N
//
//  Output:
//          state_seq       - optimal hidden state sequence given the observation
//                            sequence ob_seq
///////////////////////////////////////////////////////////////////////////////

function state_seq = hmm_viterbi(tran_prob, emit_prob, ob_seq, start_prob)
  
  N = size(tran_prob, 1);   //number of hidden states
  T = length(ob_seq);       //length of observaton sequence
    
  //delta(t,i) represents the highest probability along a single path at time t,
  //given the first t observations and hidden state in Si  
  delta = zeros(T, N);
  
  //phi is an array which keeps track of the arguments maxmizing delta for each
  //t and i  
  phi = zeros(T, N);
  
  //the hidden state sequence
  q = zeros(1, T);
  
  //initialization
  delta(1,:) = start_prob.*emit_prob(:,ob_seq(1))';
  phi(1,:) = 0;
  
  //recursion
  for t=2:T,
    for j=1:N,
      [val, idx] = max(delta(t-1,:).*tran_prob(:,j)');
      delta(t,j) = val*emit_prob(j,ob_seq(t));
      phi(t,j) = idx;
    end 
  end 
  
  //termination
  [val, idx] = max(delta(T,:));
  q(T) = idx;
    
  //path backtracking
  for t=T-1:-1:1,
    q(t) = phi(t+1, q(t+1));
  end
    
  state_seq = q;
  
endfunction