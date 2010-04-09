///////////////////////////////////////////////////////////////////////////////
// Author:  Baogang Hu <hubg@nlpr.ia.ac.cn>
// Date:    September 2009
// Version: 0.1
// Description: Calculate Normalized Mutual Information from a given m by (m+1) 
//              confusion matrix for evaluating a classifier. All NIs are 
//              calculated base on cross entropy definition. 
// Background:  Information based measures provide users for objective evaluations 
//              of classifiers. The functions below calculates NI_21 to NI_24 
//              in the references.  
// References:
//     Ref 1:   Hu, B.-G., He, R., and Yuan, X.-T., Information-Theoretic Measures 
//              for Objective Evaluation of Classifiers, submitted to a journal (2009)
//     Ref 2:   Hu, B.-G., Information Measure Toolbox for Classifier Evaluation 
//              on Open Source Software Scilab, submitted to OSSC-2009.
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
//       Input:
//             c              - Confusion matrix in size of m by (m+1),  
//                              row for exact labels,  
//                              column for prediction labels, 
//                              the (m+1)th column for rejection (or unknown) class,
//                              this matrix has to follow the constraints: 
//                                    c_ij >=0, and C_i>0 (the ith class number)
//       Output:    
//             NI             - Normalized Information listed from NI_21 to NI_24.
//                              NI_i= inf standing for singularity result 
//             A              - Accuracy. 
//             Rej            - Rejection.
//             P              - Precision for a binary classifier.
//             R              - Recall for a binary classifier.   
///////////////////////////////////////////////////////////////////////////////

function [NI,A,Rej,P,R]=confmatrix2ni_ce(c)
    ieee(2);          // = IEEE exception mode,(=0,warning and stopping when encountering singularity)
                      //   (=1, showing warning message without stopping)
                      //   (=2, without warning and stopping)
    P=[]; R=[];       // = initialization   
    n=sum(c);         // = number of total samples
    m=length(c(:,1)); // = numbers of exact classes
    Ci=sum(c,'c');    // = column vector of exact labels, m by 1
    Ci(m+1)=0;        // = adding the term for m+1 by 1 vector
    Cp=sum(c,'r');    // = row vector of prediction labels, 1 by m+1
    p_e=c/n;          // = empirical PDF of joint distribution for the two discrete variables T and Y
    pt_e=Ci'/n;       // = empirical probability mass function of T
    py_e=Cp/n;        // = empirical probability mass function of Y
    HT=0;             // = initialization of HT
    HY=0;             // = initialization of HY
    HTY=0;            // = initialization of HY of Cross Entropy
    HYT=0;            // = initialization of HY of Cross Entropy
    eps=2.2e-16;      // = error close to zero
    for i=1:m+1       
        if pt_e(i) > 0 then   // for calculation of HT         
           HT=HT-pt_e(i)*log2(pt_e(i));                               // = Eq. (3) in Ref 1
        end
        if py_e(i) > 0 then   // for calculation of HY         
           HY=HY-py_e(i)*log2(py_e(i));                               // = Eq. (3)
        end
        if py_e(i)*pt_e(i) > 0 then  // for calculation of HYT, HTY
           HTY=HTY-pt_e(i)*log2(py_e(i)); 
           HYT=HYT-py_e(i)*log2(pt_e(i));
        end
        if (py_e(i)< eps) & (pt_e(i) > 0) then
           HTY=%inf;   // label for singularity for case of ieee(0)
        end
        if (pt_e(i)< eps) & (py_e(i) > 0) then
           HYT=%inf;   // label for singularity for case of ieee(0) 
        end
    end 
    NI_21 = HT/HTY;    // Table 3 for NI_21 to NI_24
    NI_22 = HY/HYT;
    NI_23 = (HT/HTY+HY/HYT)/2;
    NI_24 = (HT+HY)/(HTY+HYT);
    NI=[NI_21 NI_22 NI_23 NI_24];
    A=sum(diag(c))/sum(c);             // Accuracy
    Rej=sum(c(:,m+1))/n                // Rejection Rate    
    if m < 3 then                      // binary classifier
         if Cp(1)>0 then
            P=c(1,1)/Cp(1);            // Precision
         else P=0;
         end
         R=c(1,1)/Ci(1);               // Recall
    end    
endfunction


//The following are one example of using the function confmatrix2ni_ce. You can remove the annotation slashes 
//before the example codes and copy the whole page into Scilab to see how function runs.
//
//    Numerical examples in the reference 
//    Examples of binary classification, Table 4
//M1=[90   0   0 ; 1   9   0];
//M2=[89   1   0 ; 0  10   0];  
//M3=[90   0   0 ; 0   9   1]; 
//M4=[89   0   1 ; 0  10   0];   
//M5=[57  38   0 ; 3   2   0];    
//M6=[89   1   0 ; 1   9   0];
//    Examples of three-class classification, Table 7
//M7 =[80  0  0  0; 0  15  0  0; 1  0  4  0 ];
//M8 =[80  0  0  0; 0  15  0  0; 0  1  4  0 ];
//M9 =[80  0  0  0; 0  15  0  0; 0  0  4  1 ];
//M10=[80  0  0  0; 1  14  0  0; 0  0  5  0 ];
//M11=[80  0  0  0; 0  14  1  0; 0  0  5  0 ];
//M12=[80  0  0  0; 0  14  0  1; 0  0  5  0 ];
//M13=[79  1  0  0; 0  15  0  0; 0  0  5  0 ];
//M14=[79  0  1  0; 0  15  0  0; 0  0  5  0 ];
//M15=[79  0  0  1; 0  15  0  0; 0  0  5  0 ];
//c=M3;
//format('v',6);
//[NI,A,Rej,P,R]=confmatrix2ni_ce(c)
