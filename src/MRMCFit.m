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
      
      vcov = table();
      reader_test_diffs = table();
      reader_means = table();

      func = @(x) summary(x, 'alpha', options.alpha);
      cells = arrayfun(func, obj.testfits, 'UniformOutput', false);
      test_means = vertcat(cells{:});

      prefix = 'summary_';
      if obj.design.reader_in_test
        prefix = strcat(prefix, 'nested_');
      end
      args = '(obj, test_means, options.alpha)';
      summary_call = @(suffix) strcat(prefix, suffix, args);

      if obj.design.fixed.id
        
        [test_equality, test_diffs] =...
          eval(summary_call('fixed_id'));
        
      elseif obj.design.fixed.reader
        
        [vcov, test_equality, test_diffs, reader_test_diffs,...
          reader_means] = eval(summary_call('fixed_reader'));
        
      else
        
        [vcov, test_equality, test_diffs] = eval(summary_call('default'));

      end
      
      res = MRMCSummary(obj, options.alpha, vcov, test_equality,...
                        test_means, test_diffs, reader_test_diffs,...
                        reader_means);
    end

    
    function comps = vcov_comps(obj, reader)

      arguments
        obj
        reader = []
      end

      tests = obj.data.test;
      readers = obj.data.reader;
      same_test = outerop(tests, tests, @eq);
      same_reader = outerop(readers, readers, @eq);

      if ~isempty(reader)
        is_group = reader == readers;
      else
        is_group = true(size(obj.data, 1), 1);
      end
      in_group = is_group * is_group';
      
      comps.n = struct('test', length(categories(obj.factors.test)),...
                       'reader', length(categories(obj.factors.reader)));
      comps.n_crosstab = size_crosstab(obj);
      comps.MS = meansq(obj);
      comps.var = mean(diag(obj.cov(is_group, is_group)));
      comps.cov = [0; 0; 0];
      f = @(x, y) mean(obj.cov(x & y & in_group));
      if obj.design.id_in_reader
        comps.cov(1) = f(~same_test, same_reader);
      elseif obj.design.id_in_test
        comps.cov(2) = f(same_test, ~same_reader);
      else
        comps.cov(1) = f(~same_test, same_reader);
        comps.cov(2) = f(same_test, ~same_reader);
        comps.cov(3) = f(~same_test, ~same_reader);
      end
    end
   
  end
  
end


function res = size_crosstab(obj)
  res = crosstab(obj.data.test, obj.data.reader);
end


function [...
  vcov, test_equality, test_diffs...
] = summary_default(obj, test_means, alpha)

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
    tinv(1 - alpha / 2, tbl.df) * tbl.StdErr * [-1 1];
  tbl.t = tbl.Estimate ./ tbl.StdErr;
  tbl{:, 'p-value'} = 2 * tcdf(abs(tbl.t), tbl.df, 'upper');
  test_diffs = tbl;

end


function [...
  vcov, test_equality, test_diffs...
] = summary_nested_default(obj, test_means, alpha)

  comps = vcov_comps(obj);
  n_test = size(comps.n_crosstab, 1);
  n_readers = sum(comps.n_crosstab, 2);
  n_reader = sum(n_readers);
  MS = comps.MS;

  n_formula = (n_reader - sum(n_readers.^2)) / (n_reader * (n_test - 1));
  denom = MS.reader + n_formula * max(comps.cov(2) - comps.cov(3), 0);

  vcov = vcov_summary(comps);

  tbl = table();
  tbl{:, 'MS(T)'} = MS.test;
  tbl{:, 'R(T)'} = MS.reader;
  tbl.Cov2 = comps.cov(2);
  tbl.Cov3 = comps.cov(3);
  tbl.Denominator = denom;
  tbl.F = MS.test / denom;
  tbl.df1 = n_test - 1;
  tbl.df2 = denom^2 / (MS.reader^2 / (n_reader - n_test));
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
  tbl.StdErr = sqrt(...
    (1 ./ n_readers(combos(:, 1)) + 1 ./ n_readers(combos(:, 2)))...
    * denom...
  );
  tbl.df = test_equality.df2;
  tbl.CI = tbl.Estimate +...
    tinv(1 - alpha / 2, tbl.df) * tbl.StdErr * [-1 1];
  tbl.t = tbl.Estimate ./ tbl.StdErr;
  tbl{:, 'p-value'} = 2 * tcdf(abs(tbl.t), tbl.df, 'upper');
  test_diffs = tbl;

end


function [...
  test_equality, test_diffs...
] = summary_fixed_id(obj, test_means, alpha)

  n = size(obj);
  MS = meansq(obj);
  
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
    tinv(1 - alpha / 2, tbl.df) * tbl.StdErr * [-1 1];
  tbl.t = tbl.Estimate ./ tbl.StdErr;
  tbl{:, 'p-value'} = 2 * tcdf(abs(tbl.t), tbl.df, 'upper');
  test_diffs = tbl;

end


function [...
  test_equality, test_diffs...
] = summary_nested_fixed_id(obj, test_means, alpha)

  sz = size_crosstab(obj);
  n_test = size(sz, 1);
  n_readers = sum(sz, 2);
  n_reader = sum(n_readers);
  MS = meansq(obj);

  tbl = table();
  tbl{:, 'MS(T)'} = MS.test;
  tbl{:, 'R(T)'} = MS.reader;
  tbl.F = MS.test / MS.reader;
  tbl.df1 = n_test - 1;
  tbl.df2 = n_reader - n_test;
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
  tbl.StdErr = sqrt(...
    (1 ./ n_readers(combos(:, 1)) + 1 ./ n_readers(combos(:, 2)))...
    * MS.reader...
  );
  tbl.df = test_equality.df2;
  tbl.CI = tbl.Estimate +...
    tinv(1 - alpha / 2, tbl.df) * tbl.StdErr * [-1 1];
  tbl.t = tbl.Estimate ./ tbl.StdErr;
  tbl{:, 'p-value'} = 2 * tcdf(abs(tbl.t), tbl.df, 'upper');
  test_diffs = tbl;

end


function [...
  vcov, test_equality, test_diffs, reader_test_diffs,...
  reader_means...
] = summary_fixed_reader(obj, test_means, alpha)

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
    norminv(1 - alpha / 2) * tbl.StdErr * [-1 1];
  tbl.z = tbl.Estimate ./ tbl.StdErr;
  tbl{:, 'p-value'} = 2 * normcdf(abs(tbl.z), 'upper');
  test_diffs = tbl;
  
  reader_test_diffs = tbl_reader_test_diffs(obj, alpha);
  
  tbl = obj.data(:, {'y', 'reader', 'test'});
  tbl.StdErr = sqrt(diag(obj.cov));
  tbl.CI = tbl.y +...
    norminv(1 - alpha / 2) * tbl.StdErr * [-1 1];
  tbl.CI = trunc(obj.y, tbl.CI);
  reader_means = tbl;

end


function [...
  vcov, test_equality, test_diffs, reader_test_diffs,...
  reader_means...
] = summary_nested_fixed_reader(obj, test_means, alpha)

  comps = vcov_comps(obj);
  n_test = size(comps.n_crosstab, 1);
  n_readers = sum(comps.n_crosstab, 2);
  n_reader = sum(n_readers);
  MS = comps.MS;

  n_formula = (n_reader - sum(n_readers.^2)) / (n_reader * (n_test - 1));
  denom = comps.var - comps.cov(2) +...
    n_formula * max(comps.cov(2) - comps.cov(3), 0);

  vcov = vcov_summary(comps);

  tbl = table();
  tbl{:, 'MS(T)'} = MS.test;
  tbl.Cov2 = comps.cov(2);
  tbl.Cov3 = comps.cov(3);
  tbl.Denominator = denom;
  tbl.X2 = (n_test - 1) * MS.test / denom;
  tbl.df = n_test - 1;
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
  tbl.StdErr = sqrt(...
    (1 ./ n_readers(combos(:, 1)) + 1 ./ n_readers(combos(:, 2)))...
    * denom...
  );
  tbl.CI = tbl.Estimate + norminv(1 - alpha / 2) * tbl.StdErr * [-1 1];
  tbl.z = tbl.Estimate ./ tbl.StdErr;
  tbl{:, 'p-value'} = 2 * normcdf(abs(tbl.z), 'upper');
  test_diffs = tbl;

  reader_test_diffs = table();

  tbl = obj.data(:, {'y', 'reader', 'test'});
  tbl.StdErr = sqrt(diag(obj.cov));
  tbl.CI = tbl.y +...
    norminv(1 - alpha / 2) * tbl.StdErr * [-1 1];
  tbl.CI = trunc(obj.y, tbl.CI);
  reader_means = tbl;

end

function res = tbl_reader_test_diffs(obj, alpha)
  n = size(obj);
  test_names = char(categories(obj.factors.test));
  reader_names = char(categories(obj.factors.reader));
  
  ests = reshape(obj.data.y, n.test, size(obj.data, 1) / n.test)';
  combos = nchoosek(1:n.test, 2);
  combos_reader = kron(ones(n.reader, 1), combos);
  reader_inds = repelem(1:n.reader, size(combos, 1));

  func = @(reader) vcov_comps(obj, reader);
  cells = arrayfun(func, reader_names);
  var = vertcat(cells.var);
  cov = horzcat(cells.cov)';
  stderrs = sqrt(2 * (var - cov(:, 1)));
  
  res = table();
  res.Reader = reader_names(reader_inds);
  res.Comparison = strcat(...
    test_names(combos_reader(:, 1)), " - ", test_names(combos_reader(:, 2))...
  );
  ests1 = ests(sub2ind(size(ests), reader_inds', combos_reader(:, 1)));
  ests2 = ests(sub2ind(size(ests), reader_inds', combos_reader(:, 2)));
  res.Estimate = ests1 - ests2;
  res.StdErr = stderrs;
  res.CI = res.Estimate + norminv(1 - alpha / 2) * res.StdErr * [-1 1];
  res.z = res.Estimate ./ res.StdErr;
  res{:, 'p-value'} = 2 * normcdf(abs(res.z), 'upper');
end


function res = tbl_test_means(obj, alpha)
  func = @(x) summary(x, 'alpha', alpha);
  cells = arrayfun(func, obj.testfits, 'UniformOutput', false);
  res = vertcat(cells{:});
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