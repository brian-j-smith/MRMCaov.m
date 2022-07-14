classdef PerformanceVariate
  
  properties (SetAccess = protected)
    label char
    truth (:, 1)
    rating (:, 1)
    keep (:, 1) logical
    metric function_handle
  end
  
  methods
    
    function obj = PerformanceVariate(truth, rating, metric)
      obj.label = 'Performance Variate';
      assert(length(truth) == length(rating),...
             'Lengths of ''truth'' and ''rating'' differ')
      obj.truth = truth;
      obj.rating = rating;
      obj.keep = true(length(truth), 1);
      obj.metric = metric;
    end
    
    
    function L = length(obj)
      L = length(obj.truth);
    end
    
    
    function plot(obj, ~)
      error(strcat("Method not defined for ", class(obj), ".")) 
    end
    
    
    function value = trunc(~, x)
      value = x;
    end
    
  end
  
end
