classdef MRMCTestFit
  
  properties
    name char;
    y PerformanceVariate;
    design struct;
    data table;
    anova struct;
    cov (:, :) double;
  end
  
  methods
    
    function obj = MRMCTestFit(name, y, design, data, anova, cov)
      obj.name = name;
      obj.y = y;
      obj.design = design;
      obj.data = data;
      obj.anova = anova;
      obj.cov = cov;
    end
    
    
    function x = mean(obj)
      x = mean(obj.data.y);
    end
    
    
    function x = meansq(obj)
      MS = obj.anova.tbl{:, 'MS'};
      x = struct('test', 0, 'reader', MS(1), 'interaction', 0);
    end
    
    function res = summary(obj, options)
      
      arguments
        obj
        options.alpha double {...
          mustBeInRange(options.alpha, 0, 1, 'exclusive')...
        } = 0.05
      end
      
      comps = vcov_comps(obj);
      n = comps.n;
      MS = comps.MS;

      res = table();
      res.Estimate = mean(obj);
      if obj.design.fixed.reader
        res{:, 'Var(Error)'} = comps.var;
        res.Cov2 = comps.cov(2);
        res.StdErr = sqrt((comps.var + (n.reader - 1) * max(res.Cov2, 0)) /...
          n.reader);
        res.CI = trunc(obj, res.Estimate +...
          res.StdErr * norminv(1 - options.alpha / 2) * [-1 1]);
      elseif obj.design.fixed.id
        res{:, 'MS(R)'} = MS.reader;
        res.StdErr = sqrt(MS.reader / n.reader);
        res.df = n.reader - 1;
        res.CI = trunc(obj, res.Estimate +...
          res.StdErr * tinv(1 - options.alpha / 2, res.df) * [-1 1]);
      else
        res{:, 'MS(R)'} = MS.reader;
        res.Cov2 = comps.cov(2);
        res.StdErr = sqrt(MS.reader / n.reader + max(res.Cov2, 0));
        res.df = (MS.reader + n.reader * max(res.Cov2, 0))^2 /...
          (MS.reader^2 / (n.reader - 1));
        res.CI = trunc(obj, res.Estimate +...
          res.StdErr * tinv(1 - options.alpha / 2, res.df) * [-1 1]);
      end
      res.Properties.RowNames = {obj.name};
            
    end
    
    
    function value = trunc(obj, x)
      value = trunc(obj.y, x);
    end
  
    
    function comps = vcov_comps(obj)
      readers = obj.data.reader;
      same_reader = outerop(readers, readers, @eq);
      
      comps.n = struct('reader', length(obj.anova.stats.n));
      comps.MS = meansq(obj);
      comps.var = mean(diag(obj.cov));
      comps.cov = [0; mean(obj.cov(~same_reader)); 0];
    end
    
  end
  
end