classdef SensitivityVariate < PerformanceVariate
  
  methods
    
    function obj = SensitivityVariate(truth, rating)
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
