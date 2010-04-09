///////////////////////////////////////////////////////////////////////////////
// Author:  Xiao-Tong Yuan
// Version: 0.1
// Date:  June 2009
// Description:  Agglomerative Mean-Shift Clustering(AggloMS)
// Reference:  Xiao-Tong Yuan, Bao-Gang Hu and Ran He, Agglomerative Mean-Shift 
//			   Clustering via Query Set Compression, SDM 2009
//
// Background: 	Mean-Shift (MS) is a powerful non-parametric clustering method. 
//				Although good accuracy can be achieved, its computational cost
//			    is particularly expensive even on moderate data sets. In this work,
//			    for the purpose of algorithm speedup, we develop an agglomerative MS
//				clustering method called Agglo-MS, along with its mode-seeking
//              ability and convergence property analysis. Our method is built
//              upon an iterative query set compression mechanism which is motivated
//              by the quadratic bounding optimization nature of MS. The whole
//              framework can be efficiently implemented in linear running time complexity.
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
//         Data           - Data matrix. Each column vector is a data point.
//         sigma          - Bandwidth of Gaussian kernel
//         ite_num        - Number of iterations.
//   Output:
//         cluster_centers        - Cluster center matrix. Each column vector is
//									a cluster center point.
//         cluster_id             - Cluster index vector.      
///////////////////////////////////////////////////////////////////////////////

function [cluster_centers, cluster_id]=aggloms(Data, sigma, ite_num)
   dim = size(Data,2);
   cluster_id = zeros(size(Data,1),1);
   Data_Query = Data;
   Data_Ref  = Data;
   
   // ----------- Kernel Density Set Compression and Clustering-------------------------
   for (i=1:ite_num)
      [sph_cov, sph_cov_id] = kerneldensitycovering(Data_Query, Data_Ref, sigma);
      if (i==1)
        cluster_id = sph_cov_id;
      else
        for (j=1:size(Data,1))
          cluster_id(j) = sph_cov_id(cluster_id(j));
        end
      end      
      Data_Query = sph_cov(:,1:dim);
   end
   cluster_centers = Data_Query;
endfunction
