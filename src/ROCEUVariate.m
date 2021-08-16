classdef ROCEUVariate < PerformanceVariate
  
  properties
    slope (1, 1) double
  end
  
  methods
    
    function obj = ROCEUVariate(truth, rating, slope)
      
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
