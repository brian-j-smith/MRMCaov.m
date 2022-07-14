classdef ROCLRposVariate < PerformanceVariate
  
  methods
    
    function obj = ROCLRposVariate(truth, rating)

% Description
%
% Designates likelihood ratio of a positive rating (LR+) as the reader
% performance metric to be computed from true case statuses and reader
% ratings for an MRMC analysis.
%
% Usage
%
%   obj = ROCLRposVariate(truth, rating)
%
% Input Arguments
%
%   truth: column array of true binary statuses.
%
%   rating: column array of numeric ratings.
%
% Details
%
%   LR+ is computed as sensitivity / (1 - specificity).
%
% Return Value
%
%   ROCLRposVariate class object.
%
% See also: mrmc

      metric = @(truth, rating) roc_lrpos(truth, rating);
      obj@PerformanceVariate(categorical(truth), categorical(rating), metric)
      assert(length(categories(obj.truth)) == 2,...
             'Argument ''truth'' must have 2 categories.')
      assert(length(categories(obj.rating)) == 2,...
             'Argument ''rating'' must have 2 categories.')
      obj.label = 'ROC Likelihood Ratio Positive Test';
    end
    
  end
  
end


function value = roc_lrpos(truth, rating)

  arguments
    truth (:, 1)
    rating(:, 1)
  end
  
  truth = categorical(truth);
  rating = categorical(rating);
  
  pos = truth == category(truth, 2);
  sens = mean(rating(pos) == category(rating, 2));
  neg = truth == category(truth, 1);
  spec = mean(rating(neg) == category(rating, 1));
  
  value = sens / (1 - spec);

end
