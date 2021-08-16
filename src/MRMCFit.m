classdef MRMCFit
  
  properties
    y PerformanceVariate
    factors table
    design struct
    data table
    anova struct
    cov (:, :) double
    testfits (1, :) MRMCTestFit
  end
  
  methods
    
    function obj = MRMCFit(y, factors, design, data, anova, cov, testfits)
      obj.y = y;
      obj.factors = factors;
      obj.design = design;
      obj.data = data;
      obj.anova = anova;
      obj.cov = cov;
      obj.testfits = testfits;
    end
    
    
    function disp(x)
      fprintf('%s ANOVA data:\n\n', class(x.y))
      disp(x.data)
      fprintf('ANOVA Table:\n\n')
      disp(x.anova.tbl(1:3, {'Source' 'd.f.' 'Sum Sq.' 'Mean Sq.'}))
      fprintf(['Obuchowski-Rockette error variance and covariance '...
               'estimates:\n\n'])
      comps = vcov_comps(x);
      tbl = table();
      tbl.Estimate = vertcat(comps.var, comps.cov);
      tbl.Correlation = tbl.Estimate / tbl.Estimate(1);
      tbl.Correlation(1) = NaN;
      tbl.Properties.RowNames = {'Error', 'Cov1', 'Cov2', 'Cov3'};
      disp(tbl)
    end
    
    
    function x = meansq(obj)
      MS = obj.anova.tbl{:, 'Mean Sq.'};
      x = struct('test', MS(2), 'reader', MS(1), 'interaction', MS(3));
    end
    
    
    function plot(obj)
      plot(obj.y, obj.factors(:, {'test', 'reader'}))
    end
    
    
    function sz = size(obj)
      varnames = {'reader' 'test'};
      sz = varfun(@(x) length(categories(x)), obj.data(:, varnames));
      sz.Properties.VariableNames = varnames;
    end
    
    
    function res = summary(obj, options)
      
      arguments
        obj
        options.alpha double {...
          mustBeInRange(options.alpha, 0, 1, 'exclusive')...
        } = 0.05
      end
      
      func = @(x) summary(x, 'alpha', options.alpha);
      cells = arrayfun(func, obj.testfits, 'UniformOutput', false);
      test_means = vertcat(cells{:});

      if obj.design.fixed.id
        
        n = size(obj);
        MS = meansq(obj);
        
        vcov = table();
        
        tbl = table();
        tbl{:, 'MS(T)'} = MS.test;
        tbl{:, 'MS(T:R)'} = MS.interaction;
        tbl.F = MS.test / MS.interaction;
        tbl.df1 = n.test - 1;
        tbl.df2 = (n.test - 1) * (n.reader - 1);
        tbl{:, 'p-value'} = fcdf(tbl.F, tbl.df1, tbl.df2, 'upper');
        test_equality = tbl;
        
        ests = test_means.Estimate;
        combos = nchoosek(1:length(ests), 2);
        test_names = categories(obj.data.test);
        tbl = table();
        tbl.Comparison = strcat(...
          test_names(combos(:, 1)), " - ", test_names(combos(:, 2))...
        );
        tbl.Estimate = ests(combos(:, 1)) - ests(combos(:, 2));
        tbl.StdErr = sqrt(2 / n.reader * MS.interaction);
        tbl.df = (n.test - 1) * (n.reader - 1);
        tbl.CI = tbl.Estimate +...
          tinv(1 - options.alpha / 2, tbl.df) * tbl.StdErr * [-1 1];
        tbl.t = tbl.Estimate / tbl.StdErr;
        tbl{:, 'p-value'} = 2 * tcdf(abs(tbl.t), tbl.df, 'upper');
        test_diffs = tbl;
        
      elseif obj.design.fixed.reader
        
        comps = vcov_comps(obj);
        n = comps.n;
        MS = comps.MS;
        denom = comps.var - comps.cov(1) + (n.reader - 1) *...
          (comps.cov(2) - comps.cov(3));
        
        vcov = vcov_summary(comps);
        
        tbl = table();
        tbl{:, 'MS(T)'} = MS.test;
        tbl.Cov1 = comps.cov(1);
        tbl.Cov2 = comps.cov(2);
        tbl.Cov3 = comps.cov(3);
        tbl.Denominator = denom;
        tbl.X2 = (n.test - 1) * MS.test / denom;
        tbl.df = n.test - 1;
        tbl{:, 'p-value'} = chi2cdf(tbl.X2, tbl.df, 'upper');
        test_equality = tbl;

        ests = test_means.Estimate;
        combos = nchoosek(1:length(ests), 2);
        test_names = categories(obj.data.test);
        tbl = table();
        tbl.Comparison = strcat(...
          test_names(combos(:, 1)), " - ", test_names(combos(:, 2))...
        );
        tbl.Estimate = ests(combos(:, 1)) - ests(combos(:, 2));
        tbl.StdErr = sqrt(2 / n.reader * denom);
        tbl.CI = tbl.Estimate +...
          norminv(1 - options.alpha / 2) * tbl.StdErr * [-1 1];
        tbl.t = tbl.Estimate / tbl.StdErr;
        tbl{:, 'p-value'} = 2 * normcdf(abs(tbl.t), 'upper');
        test_diffs = tbl;

      else
        
        comps = vcov_comps(obj);
        n = comps.n;
        MS = comps.MS;
        denom = MS.interaction + n.reader * max(comps.cov(2) - comps.cov(3), 0);

        vcov = vcov_summary(comps);

        tbl = table();
        tbl{:, 'MS(T)'} = MS.test;
        tbl{:, 'MS(T:R)'} = MS.interaction;
        tbl.Cov2 = comps.cov(2);
        tbl.Cov3 = comps.cov(3);
        tbl.Denominator = denom;
        tbl.F = MS.test / denom;
        tbl.df1 = n.test - 1;
        tbl.df2 = denom^2 /...
          (MS.interaction^2 / ((n.test - 1) * (n.reader - 1)));
        tbl{:, 'p-value'} = fcdf(tbl.F, tbl.df1, tbl.df2, 'upper');
        test_equality = tbl;

        ests = test_means.Estimate;
        combos = nchoosek(1:length(ests), 2);
        test_names = categories(obj.data.test);
        tbl = table();
        tbl.Comparison = strcat(...
          test_names(combos(:, 1)), " - ", test_names(combos(:, 2))...
        );
        tbl.Estimate = ests(combos(:, 1)) - ests(combos(:, 2));
        tbl.StdErr = sqrt(2 / n.reader * denom);
        tbl.df = test_equality.df2;
        tbl.CI = tbl.Estimate +...
          tinv(1 - options.alpha / 2, tbl.df) * tbl.StdErr * [-1 1];
        tbl.t = tbl.Estimate / tbl.StdErr;
        tbl{:, 'p-value'} = 2 * tcdf(abs(tbl.t), tbl.df, 'upper');
        test_diffs = tbl;
        
      end
      
      res = MRMCSummary(obj, options.alpha, vcov, test_equality, test_means,...
                        test_diffs);
    end

    
    function comps = vcov_comps(obj)
      tests = obj.data.test;
      readers = obj.data.reader;
      same_test = outerop(tests, tests, @eq);
      same_reader = outerop(readers, readers, @eq);
      
      comps.n = struct('test', length(categories(obj.factors.test)),...
                       'reader', length(categories(obj.factors.reader)));
      comps.MS = meansq(obj);
      comps.var = mean(diag(obj.cov));
      comps.cov = [mean(obj.cov(~same_test & same_reader));...
                   mean(obj.cov(same_test & ~same_reader));...
                   mean(obj.cov(~same_test & ~same_reader))];
    end
   
  end
  
end


function res = vcov_summary(comps)

  MS = comps.MS;
  n = comps.n;
  
  res = table();
  res.Estimate = vertcat(...
    (MS.reader - MS.interaction) / n.test - comps.cov(1) +...
      comps.cov(3),...
    MS.interaction - comps.var + comps.cov(1) +...
      (comps.cov(2) - comps.cov(3)),...
    comps.var, comps.cov...
  );
  res.Correlation = res.Estimate / comps.var;
  res.Correlation(1:3) = NaN;
  res.Properties.RowNames = {'reader' 'reader*test' 'Error' 'Cov1'...
                              'Cov2' 'Cov3'};

end