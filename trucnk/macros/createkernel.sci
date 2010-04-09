///////////////////////////////////////////////////////////////////////////////
// Author:  Jia Wu
// Version: 0.1
// Date: Nov 2009
// Description: Kernel Function
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
//        x       - dim*numx data matrix. Each column is a data point.
//        y       - dim*numy data matrix. Each column is a data point.
//        param   - A struct variable with the following fields:
//                  typ  - Gaussian:      exp(-|x-y|^2/2t^2)
//                         Polynomial:    (c*x'*y+r)^d
//                         Linear:         x'*y
//                         Sigmoid:        tanh(c*x'*y+r)
//                  t    - kernel parameter
//                  c    - kernel parameter
//                  r    - kernel parameter
//                  d    - kernel parameter
//		
//  Output:
//        K       - numx*numy kernel matrix. 
///////////////////////////////////////////////////////////////////////////////

function K = createkernel(x, y, param)
  
  if ~isstruct(param),
    error('param should be a struct variable.');
  end
  
  if isempty(y),
    y = x;
  end
    
  select param.typ
    
    case 'Gaussian' then
      if isempty(param.t),
        param.t = 1;
      end
      
      r = size(x, 2);
      c = size(y, 2);
      K = zeros(r, c);
      for i = 1:r,
        for j = 1:c,
          K(i, j) = norm(x(:,i)-y(:,j));
        end
      end
      K = exp(-K.^2/(2*param.t^2));
          
    case 'Polynomial' then
      if isempty(param.c)
        param.c = 1;
      end
      if isempty(param.r)
        param.r = 0;
      end
      if isempty(param.d)
        param.d = 2;
      end
      
      K = (param.c*(x'*y)+param.r).^param.d;
  
    case 'Linear' then
      K = x'*y;
    
    case 'Sigmoid' then
      if isempty(param.c)
        param.c = 1;
      end
      if isempty(param.r)
        param.r = 0;
      end
      
      K = (param.c*x'*y+param.r);
    
    else
      K = [];
      error('No kernel type assigned, or wrong kernel type.');
  end

endfunction
