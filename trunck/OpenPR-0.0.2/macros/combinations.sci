///////////////////////////////////////////////////////////////////////////////
// Author:   Jia Wu
// Date:     June 2010
// Description:  find all the combinations of given number of indices 
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
//         indices      - input indices and should be a vector
//         number       - the required number of indices
//
//  Output:
//         comb_mat     - all the combinations of the required number of indices
//                        the size is comb_num*number; comb_num is the numbers of
//                        all possible combinations, and num is the required number
//                        of indices 
///////////////////////////////////////////////////////////////////////////////

function comb_mat = combinations(indices, number)
  
  indices = indices(:);
  
  idx_num = length(indices);
  
  //combination_num = factorial(idx_num)/(factorial(number)*factorial(idx_num-number));
  
  if number==idx_num,
    comb_mat = indices';
  elseif number==1,
    comb_mat = indices;
  else
    //recursive
    comb_mat = [];
    for i=1:idx_num-number+1,
      temp = combinations(indices(i+1:idx_num), number-1);
      comb_mat = [comb_mat; indices(i)*ones(size(temp, 1), 1) temp];
    end
  end
   
endfunction
