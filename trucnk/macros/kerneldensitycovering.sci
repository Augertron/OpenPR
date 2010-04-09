///////////////////////////////////////////////////////////////////////////////
// Author:  Xiao-Tong Yuan
// Version: 0.1
// Date:  June 2009
// Description: This function is used to construct a set of hyperspheres to cover
//			    the query dataset.
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
//         Data_Query     - Query data matrix. Each column vector is a data point.
//         Data_Ref       - Reference data matrix. Each column vector is a data point.
//         sigma          - Bandwidth of Gaussian kernel
//
//   Output:
//         sph_cov        - Constructed hyperspheres that cover the query data set.
//         sph_cov_id     - Vector of covering hypersphere index. Each data point
//							receives one index.   
///////////////////////////////////////////////////////////////////////////////

function [sph_cov, sph_cov_id]=kerneldensitycovering(Data_Query, Data_Ref, sigma)
  dim = size(Data,2);
  sph_cov = [];
  sph_cov_count = 0;
  sph_cov_id = [];
  
  //------------- Scan the data set and sequentially cover it with hyperspheres ----------------
  for (i=1:size(Data_Query,1))
       data_cur = Data_Query(i,:);
       cover_flag = 0;
       
       for (j=1:size(sph_cov,1)) // check whether the current point is covered by any exisitng hyperspheres
        sph_cur_center = sph_cov(j,1:dim);
        sph_cur_radius = sph_cov(j,dim+1);
        dist = sqrt(sum((data_cur-sph_cur_center).^2,2));
         if ((dist<sph_cur_radius)|(dist<1e-3))
              sph_cov_id(i) = j; 
              cover_flag = 1;
              break;
         end
       end    
       if (cover_flag==0) // if the current has not been covered yet, construct a new hypersphere
            [s_center, s_radius] = sphereconstruction(data_cur, Data_Ref, sigma); 
            sph_cov = [sph_cov; [s_center, s_radius]];
            sph_cov_count = sph_cov_count+1;
            sph_cov_id(i) = sph_cov_count;
       end
              
end
endfunction
