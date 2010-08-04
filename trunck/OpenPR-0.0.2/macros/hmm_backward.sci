///////////////////////////////////////////////////////////////////////////////
// Author:   Jia Wu
// Date:     July 2010
// Description:  backward algorithm for computing the probabilty of an observation
//               sequence given a hidden markov model 
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
//          P               - probability of an observation sequence generated
//                            by the given model 
//          be              - probability of matrix of size T*N
//                            be(t,i) represents the probability of remainder
//                            partial observation sequence O(t+1)O(t+1)...O(T) 
//                            given hidden state Si at time t
///////////////////////////////////////////////////////////////////////////////

function [P, be] = hmm_backward(tran_prob, emit_prob, ob_seq, start_prob)
  
  N = size(tran_prob, 1);   //number of hidden states
  M = size(emit_prob, 2);   //number of visible states
  T = length(ob_seq);       //length of observaton sequence
  
  //backward variables
  be = zeros(T, N);
  
  //initialization
  be(T,:) = 1;
  
  //induction
  for t=T-1:-1:1,
    for i=1:N,
      be(t,i) = sum(tran_prob(i,:).*emit_prob(:,ob_seq(t+1))'.*be(t+1,:));
    end
  end
    
  //termination
  P = sum(start_prob.*emit_prob(:,ob_seq(1))'.*be(1,:));  
  
endfunction

