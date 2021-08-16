classdef ROCLRposVariate < PerformanceVariate
  
  methods
    
    function obj = ROCLRposVariate(truth, rating)
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
