///////////////////////////////////////////////////////////////////////////////
// Author:   Jia Wu
// Date:     July 2010
// Description:  foward algorithm for computing the probabilty of an observation
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
//          alpha           - probability matrix of size T*N 
//                            alpha(t,i) represents the probability of partial
//                            observation sequence O(1)O(2)...O(t) and hidden
//                            state being Si at time t 
///////////////////////////////////////////////////////////////////////////////

function [P, alpha] = hmm_forward(tran_prob, emit_prob, ob_seq, start_prob)
  
  N = size(tran_prob, 1);   //number of hidden states
  M = size(emit_prob, 2);   //number of visible states
  T = length(ob_seq);       //length of observaton sequence
  
  //forward variables(probability of partial observation sequence)
  alpha = zeros(T, N);
  
  //initialization
  alpha(1,:) = start_prob.*emit_prob(:,ob_seq(1))';
  
  //induction
  for t=2:T,
    for j=1:N,
      alpha(t,j) = (alpha(t-1,:)*tran_prob(:,j))*emit_prob(j,ob_seq(t));
    end
  end     
  
  //termination
  P = sum(alpha(T,:));
  
endfunction