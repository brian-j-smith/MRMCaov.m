classdef ROCEUVariate < PerformanceVariate
  
  properties
    slope (1, 1) double
  end
  
  methods
    
    function obj = ROCEUVariate(truth, rating, slope)
 
% Description
%
% Designates expected utility of the ROC curve as the reader performance
% metric to be computed from true case statuses and reader ratings for an
% MRMC analysis.
%
% Usage
%
%   obj = ROCEUVariate(truth, rating, slope)
%
% Input Arguments
%
%   truth: column array of true binary statuses.
%
%   rating: column array of numeric ratings.
%
%   slope: numeric slope value at which to compute expected utility
%     [default: 1].
%
% Return Value
%
%   ROCEUVariate class object.
%
% References
%
% Abbey CK, Samuelson FW and Gallas BD (2013). Statistical power
% considerations for a utility endpoint in observer performance studies.
% Academic Radiology, 20: 798-806.
%
% See also: mrmc
      
      arguments
        truth
        rating
        slope = 1
      end
      
      metric = @(truth, rating) roc_eu(truth, rating, slope);
      obj@PerformanceVariate(categorical(truth), double(rating), metric)
      assert(length(categories(obj.truth)) == 2,...
             'Argument ''truth'' must have 2 categories.')
      obj.label = 'ROC Expected Utility';
      obj.slope = slope;
    end
    
  end
  
end


function value = roc_eu(truth, rating, slope)

  arguments
    truth (:, 1)
    rating(:, 1) double
    slope (1, 1) double = 1
  end
  
  truth = categorical(truth);
  [fpr, tpr] = perfcurve(truth, rating, category(truth, 2));
  value = max(tpr - slope * fpr);

end
