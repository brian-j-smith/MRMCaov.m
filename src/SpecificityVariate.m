classdef SpecificityVariate < PerformanceVariate
  
  methods
    
    function obj = SpecificityVariate(truth, rating)
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
