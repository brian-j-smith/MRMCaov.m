classdef SensitivityVariate < PerformanceVariate
  
  methods
    
    function obj = SensitivityVariate(truth, rating)

% Description
%
% Designates sensitivity as the reader performance metric to be computed
% from true case statuses and reader ratings for an MRMC analysis.
%
% Usage
%
%   obj = SensitivityVariate(truth, rating)
%
% Input Arguments
%
%   truth: column array of true binary statuses.
%
%   rating: column array of binary ratings.
%
% Details
%
%   Sensitivity is calculated from the data as the proportion of positive
%   cases who have positive ratings.  As such it will not be estimable for
%   reader-test combinations that do not have any positive cases.
%
% Return Value
%
%   SensitivityVariate class object.
%
% See also: mrmc

      metric = @(truth, rating) mean(rating == category(rating, 2));
      obj@PerformanceVariate(categorical(truth), categorical(rating), metric)
      assert(length(categories(obj.truth)) == 2,...
             'Argument ''truth'' must have 2 categories.')
      assert(length(categories(obj.rating)) == 2,...
             'Argument ''rating'' must have 2 categories.')
      obj.label = 'Sensitivity';
      obj.keep = obj.truth == category(obj.truth, 2);
      obj.truth = obj.truth(obj.keep);
      obj.rating = obj.rating(obj.keep);
    end
    
    function value = trunc(~, x)
      value = max(0, min(x, 1));
    end
    
  end
  
end
