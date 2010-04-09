///////////////////////////////////////////////////////////////////////////////
// Author:  Xiao-Tong Yuan
// Version: 0.1
// Date:  June 2009
// Description: This function is used to construct a hypersphere from a given 
//				point by running one-iteration of Mean-Shift on the data set.
// Reference:   Xiao-Tong Yuan, Bao-Gang Hu and Ran He, Agglomerative Mean-Shift
//				Clustering via Query Set Compression, SDM 2009
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
//   Input:
//         p_start        - Given point from which hypersphere to be constructed
//         data           - Data matrix. Each column vector is a data point.
//         sigma          - Bandwidth of Gaussian kernel
//
//   Output:
//         s_center       - Center of the constructed hypersphere.
//         radius         - Radius of the cosntructed hypersphere.  
///////////////////////////////////////////////////////////////////////////////

function [s_center, radius]=sphereconstruction(p_start, data, sigma)
 // For super sphere construction
  S_n = 0;
  S_d = 0;
  
  p_start_rep = ones(size(data,1),1)*p_start;
 
  weight = exp(-sum((p_start_rep-data).^2,2)/sigma^2);
  weight_rep = weight*ones(1, size(data,2));
  S_n = sum(weight_rep.*data,1);
  S_d = sum(weight,1);


s_center = S_n/S_d;
radius = sqrt(sum((s_center-p_start).^2,2));
endfunction
