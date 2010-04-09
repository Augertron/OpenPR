///////////////////////////////////////////////////////////////////////////////
// Author:  Ran He
// Version: 0.1
// Date: May 2009
// Description: Weighted principal component analysis(WPCA)
// Reference: C.M. Bishop. Pattern Recognition and Machine Learning. Information Science and Statistics,2006
//
// Background: Principal Component Analysis (PCA) is a linear data transformation technique which 
//             plays an important role in the studies of computer vision and machine learning. It 
//             has been widely used for the representation of high dimensional data such as appea
//             -rance, shape, visual tracking, etc. And it is commonly used as a preprocessing st
//             -ep to project high-dimensional data into a low-dimensional subspace and reduce the 
//             noise at same time.
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
//        Train_Patterns - Data matrix. Each column vector is a data point.
//        weight         - A column vector whose length is equal to number 
//                         of data points.
//  Output:
//        eig_vec        - Each column is an eigvector. eig_vec'*eig_vec=I.
//        m_vec          - data center. 
//        eig_val        - The sorted eigvalue of WPCA algorithm.
///////////////////////////////////////////////////////////////////////////////

function [eig_vec,m_vec,eig_val] = wpca(Train_Patterns,weight)
      
    [dim, sample_num] = size(Train_Patterns);
    m_vec = Train_Patterns*(weight/sum(weight));
    for i=1:sample_num
    	Train_Patterns(:,i)= sqrt(weight(i))*(Train_Patterns(:,i)- m_vec);
    end
    
    if (dim<=sample_num) then
        cov_matrix = (1./sum(weight))*Train_Patterns*Train_Patterns';
        [eig_vec, eig_val] = spec(cov_matrix);   
        
        //sort the eigvector    
        [eig_val, I] = sort(diag(eig_val));
        eig_vec = eig_vec(:,I);
    
    else // see Bishop' book for details
        cov_matrix = (1./sum(weight))*Train_Patterns'*Train_Patterns;
        [eig_vec, eig_val] = spec(cov_matrix);
        eig_val = diag(eig_val);
        [eig_val, I] = sort(eig_val);
        eig_vec = eig_vec(:,I);
        eig_val = eig_val(1:$-1);
        eig_vec = eig_vec(:,1:$-1);
        eig_vec = Train_Patterns*eig_vec*(sum(weight).^(-0.5))*diag(eig_val.^(-0.5)); 
    end

    //normalize the eigvector
    for i = 1:size(eig_vec,2)
    	eig_vec(:,i) = eig_vec(:,i)./norm(eig_vec(:,i));
    end
endfunction

//The following are two examples of using the function wpca. You can remove the annotation slashes 
//before the example codes and copy the whole page into Scilab to see how to use the wpca function.  
//
//e.g.1: learn the principal component of a dataset
//
//         Train_Patterns = rand(2,100);   //number of dimension is 2 number of data points is 100
//         weight = ones(100,1);           //A column vector of weight
//         [eig_vec,m_vec,eig_val] = wpca(Train_Patterns,weight);
//
//e.g.2: a graph example of principal component
//          
//		   rotation = [7 -cos(3.14/4);sin(3.14/4) 1]
//         Train_Patterns = rand(2,100);
//         Train_Patterns = rotation* Train_Patterns;                  
//         weight = ones(100,1);           //equal weight
//         plot2d(Train_Patterns(1,:),Train_Patterns(2,:),style=-4); //plot the data points
//         [eig_vec,m_vec,eig_val] = wpca(Train_Patterns,weight);
//         plot2d(m_vec(1),m_vec(2),style=-5);                       //plot the mean vector
//         //plot the first principal component
//         x=-1:0.1:7;
//         y=(eig_vec(2,1)/eig_vec(1,1))*(x-m_vec(1))+m_vec(2);
//         plot2d(x,y,style=2);
//         legends(["data";"data center";"eigen vector"],[-4,-5,2], with_box=%f, opt="?")
