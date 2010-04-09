///////////////////////////////////////////////////////////////////////////////
// Author:  Baogang Hu <hubg@nlpr.ia.ac.cn>
// Date:    September 2009
// Version: 0.1
// Description: Calculate Normalized Mutual Information from a given m by (m+1)
//              confusion matrix for evaluating a classifier. All NIs are
//              calculated base on mutual information definition. 
// Background:  Information based measures provide users for objective evaluations
//              of classifiers. The function calculates NI_1 to NI_9 in the references.                       
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
//       Output:    NI,A,Rej,P,R
//             NI             - Normalized Information from NI_1 to NI_9.
//             A              - Accuracy. 
//             Rej            - Rejection.
//             P              - Precision for a binary classifier.
//             R              - Recall for a binary classifier.
///////////////////////////////////////////////////////////////////////////////

function [NI,A,Rej,P,R]=confmatrix2ni_mi(c)
    ieee(2);          // = IEEE exception mode,(=0,warning and stopping when encountering singularity)
                      //   (=1, showing warning message without stopping)
                      //   (=2, without warning and stopping)
    P=[]; R=[];       // = initialization   
    n=sum(c);         // = number of total samples
    m=length(c(:,1)); // = number of exact classes
    m=int32(m);       // = transfer into integer 
    Ci=sum(c,'c');    // = column vector of exact labels, each term has to be > 0
    Cp=sum(c,'r');    // = row vector of prediction labels, 1 by m+1
    p_e=c/n;          // = empirical PDF of joint distribution for the two discrete variables T and Y
    pt_e=Ci'/n;       // = empirical probability mass function of T
    py_e=Cp/n;        // = empirical probability mass function of Y
    MI=0;             // = initialization of MI
    HT=0;             // = initialization of HT
    HY=0;             // = initialization of HY
    eps=2.2e-16;      // = error close to zero
    for j=1:m+1       // for calculation of MI
        for i=1:m
            if p_e(i,j) > 0 then
               if pt_e(i)*py_e(j) > 0 then            
                  MI=MI+p_e(i,j)*log2(p_e(i,j)/pt_e(i)/(py_e(j)));  // = Eq. (4) in Ref 1
               end
            end
            if p_e(i,j) > 0 then
               
            end
            if (i==m) & (j==m) then MI_m=MI; end                    // = summation of j up to m instead of m+1, modified MI on Eq. (4)
        end
    end
    for i=1:m         // for calculation of HT
        HT=HT-pt_e(i)*log2(pt_e(i));                                // = Eq. (3)
    end
    for j=1:m+1       // for calculation of HY
        if py_e(j) > 0 then            
           HY=HY-py_e(j)*log2(py_e(j));                             // = Eq. (3)
        end
    end
    MI=clean(MI);     // round to zero for very small error entries
    MI_m=clean(MI_m);   
    NI_1 = MI/HT;     // Table 1 for NI_1 to NI_9 
    NI_2 = MI_m/HT;
    if HY < eps then 
       NI_3 = %inf;   // = singular for case of ieee(0)
       NI_4 = %inf;
       NI_5 = %inf;
       NI_6 = %inf;
       NI_9 = %inf;
    else  
       NI_3 = MI/HY;
       NI_4 = (MI/HT+MI/HY)/2;
       NI_5 = 2*MI/(HT+HY);
       NI_6 = MI/sqrt(HT*HY);
       NI_9 = MI/min([HT,HY]);
    end
    NI_7 = MI/(HT+HY-MI);
    NI_8 = MI/max([HT,HY]);
    NI=[NI_1 NI_2 NI_3 NI_4 NI_5 NI_6 NI_7 NI_8 NI_9];
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


//The following are one example of using the function confmatrix2ni_mi. You can remove the annotation slashes 
//before the example codes and copy the whole page into Scilab to see how function runs.
//  
//    Numerical examples in the reference
//    Examples of binary classification, Table 4 in Ref 1
//M1=[90   0   0 ; 1   9   0];
//M2=[89   1   0 ; 0  10   0];  
//M3=[90   0   0 ; 0   9   1]; 
//M4=[89   0   1 ; 0  10   0];   
//M5=[57  38   0 ; 3   2   0];    
//M6=[89   1   0 ; 1   9   0];
//    Examples of three-class classification, Table 7 in Ref 1
//M7 =[80  0  0  0; 0  15  0  0; 1  0  4  0 ];
//M8 =[80  0  0  0; 0  15  0  0; 0  1  4  0 ];
//M9 =[80  0  0  0; 0  15  0  0; 0  0  4  1 ];
//M10=[80  0  0  0; 1  14  0  0; 0  0  5  0 ];
//M11=[80  0  0  0; 0  14  1  0; 0  0  5  0 ];
//M12=[80  0  0  0; 0  14  0  1; 0  0  5  0 ];
//M13=[79  1  0  0; 0  15  0  0; 0  0  5  0 ];
//M14=[79  0  1  0; 0  15  0  0; 0  0  5  0 ];
//M15=[79  0  0  1; 0  15  0  0; 0  0  5  0 ];
//c=M1;
//format('v',6);
//[NI,A,Rej,P,R]=confmatrix2ni_mi(c)
