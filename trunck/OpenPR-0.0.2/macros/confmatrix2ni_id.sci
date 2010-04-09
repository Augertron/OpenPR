///////////////////////////////////////////////////////////////////////////////
// Author:  Baogang Hu <hubg@nlpr.ia.ac.cn>
// Date:    September 2009
// Version: 0.1
// Description: Calculate Normalized Mutual Information from a given m by (m+1)
//              confusion matrix for evaluating a classifier. All NIs are
//              calculated base on information divergence definition. 
// Background:  Information based measures provide users for objective evaluations 
//              of classifiers. The function below calculates NI_10 to NI_20 
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
//             NI             - Normalized Information listed from NI_10 to NI_20.
//                              NI_i= inf standing for singularity result 
//             A              - Accuracy. 
//             Rej            - Rejection.
//             P              - Precision for a binary classifier.
//             R              - Recall for a binary classifier.   
//
///////////////////////////////////////////////////////////////////////////////

function [NI,A,Rej,P,R]=confmatrix2ni_id(c)
    ieee(2);          // = IEEE exception mode,(=0,warning and stopping when encountering singularity)
                      //   (=1, showing warning message without stopping)
                      //   (=2, without warning and stopping)
    P=[]; R=[];       // = initialization
    n=sum(c);         // = number of total samples
    m=length(c(:,1)); // = numbers of exact classes
    Ci=sum(c,'c');    // = column vector of exact labels, m by 1
    Ci(m+1)=0;        // = adding the term for m+1 by 1 vector
    Cp=sum(c,'r');    // = row vector of prediction labels, 1 by m+1
    p=Ci'/n;          // = empirical probability mass function of T
    q=Cp/n;           // = empirical probability mass function of Y
    r=(p+q)/2;        // = means 
    eps=2.2e-16;      // = error close to zero
    QMI=0;            // = initialization
    CS1=0;
    CS2=0;
    CS3=0;
    KL=0;             // = initialization 
    KLI=0;            // 
    JSpr=0;           // = initialization  
    JSqr=0;           // = initialization  
    J=0;              // = initialization  
    KLQ=0;            // = initialization  
    VD=0;             // = initialization  
    HD=0;             // = initialization  
    BD=0;             // = initialization  
    X2p=0;            // = initialization  
    X2q=0;            // = initialization  
    IS1=0;            // = initialization  
    IS2=0;            // = initialization  
    for i=1:m+1
      QMI=QMI+(p(i)-q(i))^2;
      CS1=CS1+p(i)^2;
      CS2=CS2+q(i)^2;
      CS3=CS3+p(i)*q(i);
      if p(i)> 0 then
        if q(i)> 0 then
           KLpq=p(i)*log2(p(i)/q(i));
        else
           KLpq=-%inf;  // for case of ieee(0)
        end
        KLpr=p(i)*log2(p(i)/r(i));
        X2p= X2p+(p(i)-q(i))^2/p(i);
      else
        KLpq=0;
        KLpr=0;
        if q(i) > 0 then
           X2p=%inf;    // for case of ieee(0)
        end
      end
      if q(i)> 0 then
        if p(i)> 0 then
           KLqp=q(i)*log2(q(i)/p(i));   
        else
           KLqp=-%inf;  // for case of ieee(0)
        end
        KLqr=q(i)*log2(q(i)/r(i));
        X2q= X2q+(p(i)-q(i))^2/q(i);
      else
        KLqp=0;
        KLqr=0;
        if p(i) > 0 then
           X2q=%inf;    // for case of ieee(0)
        end
      end
      KLI=KLI+KLqp;
      KL=KL+KLpq;
      JSpr=JSpr+KLpr;
      JSqr=JSqr+KLqr;
      KLQ=KLQ+(p(i)-q(i))^2;
      VD=VD+abs(p(i)-q(i));
      HD=HD+(sqrt(p(i))-sqrt(q(i)))^2;
      BD= BD+sqrt(p(i)*q(i));
    end
    CS=log2(CS1*CS2/CS3^2);
    J=KLI+KL;
    NI_11=exp(-CS);        // CS-Quadratic Divergence
    if abs(J)<eps then     // Singularity checking
      RA=%inf;             // for case of ieee(0)
    else
      RA=KLI*KL/J;      
    end
    JS=(JSpr+JSqr);
    BD=-log2(BD);
    SX2=X2p+X2q;
    NI_10=exp(-QMI);       // ED-Quadratic Divergence 
    NI_11=exp(-CS);        // CS-Quadratic Divergence
    if abs(KL)==%inf then  // Singularity checking
      NI_12=%inf;          // for case of ieee(0)
    else
      NI_12=exp(-KL);      // KL Divergence
    end
    NI_13=exp(-BD);        // Bhattacharyya Distance
    if abs(X2q)==%inf then // Singularity checking
      NI_14=%inf;          // for case of ieee(0)
    else
      NI_14=exp(-X2q);     // X2 (Pearson) Divergence
    end
        NI_15=exp(-HD);    // Hellinger Distance
    NI_16=exp(-VD);        // Variation Distance
    if abs(J)==%inf then   // Singularity checking
      NI_17=%inf;          // for case of ieee(0)
    else
      NI_17=exp(-J);       // J divergence (Symmetric KL divergence)
    end
    if abs(JS)==%inf then  // Singularity checking
      NI_18=%inf;          // for case of ieee(0)
    else
      NI_18=exp(-JS);      // L (or JS) divergence
    end
    if abs(SX2)==%inf then // Singularity checking
      NI_19=%inf;          // for case of ieee(0)
    else
      NI_19=exp(-SX2);     // Symmetric X2 Divergence
    end
    if ((abs(RA)==%inf) | (string(RA)=='Nan')) then 
      NI_20=%inf;
    else
      NI_20=exp(-RA);      // Resistor Average Distance
    end 
    NI=[NI_10 NI_11 NI_12 NI_13 NI_14 NI_15 NI_16 NI_17 NI_18 NI_19 NI_20];
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


//The following are one example of using the function confmatrix2ni_id. You can remove the annotation slashes 
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
//c=M1;
//format('v',7);
//[NI,A,Rej,P,R]=confmatrix2ni_id(c)
