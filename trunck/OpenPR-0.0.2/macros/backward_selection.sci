///////////////////////////////////////////////////////////////////////////////
// Author:   Jia Wu
// Date:     June 2010
// Description:  backward feature selection
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
//          samples       - dim*num data matrix;
//                          each column is a data point, each row is a feature
//          labels        - a 1*num vector of labels for each data point
//          final_dim     - the number of output dimension after feature selection
//          fold          - number of folds for cross validation
//
//  Output:
//          new_samples   - final_dim*num data matrix after feature selection
//          feature_idx   - selected feature indices
///////////////////////////////////////////////////////////////////////////////

function [new_samples, feature_idx] = backward_selection(samples, labels, final_dim, fold)
  
  [dim, num] = size(samples);
  
  if(final_dim > dim)
    error('The number of feature dimension after feature selection should be smaller than the original dimension.');
  end
  
  //construct cross validation folds
  tmp_number = floor(num/fold)*fold; 
  tmp = matrix([1:tmp_number], fold, tmp_number/fold);
  train_idx = zeros(fold, tmp_number/fold*(fold-1));
  test_idx = zeros(fold, tmp_number/fold);
  for i=1:fold
    train_idx(i, :) = matrix(tmp([1:i-1,i+1:fold], :), 1, tmp_number/fold*(fold-1));
    test_idx(i, :) = tmp(i, :);
  end
  
  //initial feature subset
  subsets = 1:dim;
  
  for i=dim:-1:final_dim,
    com_num = size(subsets, 1);
    score = zeros(fold, com_num);
    for j=1:com_num,
      for k=1:fold,
        train_samples = samples(subsets(j,:), train_idx(k,:));
        train_labels = labels(train_idx(k,:));
        test_samples = samples(subsets(j,:), test_idx(k,:));
        
        //train the classifier using the train_samples        
        //classify test_samples
        //using least-squares algorithm (classifiers could be chosen later) 
        pre_labels = least_squares(train_samples, train_labels, test_samples);
        
        //calculate score using the test_samples
        score(k, j) = 1-mean(bitxor(pre_labels, labels(test_idx(k, :))));
      end
    end
    
    //find the best feature subset
    avg_score = mean(score, 'r');
    [val, best_subset] = max(avg_score);
    select_feature = subsets(best_subset, :);   
    
    if i==final_dim,
      break;
    else    
      //delete one feature from the select_feature
      comb_mat = combinations(1:length(select_feature), length(select_feature)-1);
      //new feature subsets
      subsets = zeros(comb_mat);
      for ii=1:size(comb_mat,1),
        subsets(ii, :) = select_feature(comb_mat(ii,:));
      end  
    end       
    
  end   
  
  //new samples points composed of the selected features
  new_samples = samples(select_feature, :);
  //selected features
  feature_idx = select_feature;   

endfunction