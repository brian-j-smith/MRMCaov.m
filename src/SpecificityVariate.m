classdef SpecificityVariate < PerformanceVariate
  
  methods
    
    function obj = SpecificityVariate(truth, rating)

% Description
%
% Designates specificity as the reader performance metric to be computed
% from true case statuses and reader ratings for an MRMC analysis.
%
% Usage
%
%   obj = obj = SpecificityVariate(truth, rating)
%
% Input Arguments
%
%   truth: column array of true binary statuses.
%
%   rating: column array of binary ratings.
%
% Details
%
%   Specificity is calculated from the data as the proportion of negative
%   cases who have negative ratings.  As such it will not be estimable for
%   reader-test combinations that do not have any negative cases.
%
% Return Value
%
%   SpecificityVariate class object.
%
% See also: mrmc

      metric = @(truth, rating) mean(rating == category(rating, 1));
      obj@PerformanceVariate(categorical(truth), categorical(rating), metric)
      assert(length(categories(obj.truth)) == 2,...
             'Argument ''truth'' must have 2 categories.')
      assert(length(categories(obj.rating)) == 2,...
             'Argument ''rating'' must have 2 categories.')
      obj.label = 'Specificity';
      obj.keep = obj.truth == category(obj.truth, 1);
      obj.truth = obj.truth(obj.keep);
      obj.rating = obj.rating(obj.keep);
    end
    
    function value = trunc(~, x)
      value = max(0, min(x, 1));
    end
    
  end
  
end
