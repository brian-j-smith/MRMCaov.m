function fit = mrmc(y, test, reader, id, options)

% Multi-Reader Multi-Case Statistical Analysis
% 
% Description
%
%   Computes performance metrics and covariances needed to perform a
%   multi-reader multi-case statistical analyis.
%
% Usage
%
%   fit = mrmc(y, test, reader, id, options)
%
% Input Arguments
%
%   y: PerformanceVariate object defining true case statuses, corresponding
%     reader ratings, and a reader performance metric to compute on them.
%
%   test, reader, id: column arrays of grouping variables that identify the
%     test modality, reader, and case for the observations in y.  Each of
%     the variables can be either categorical, numerical, character, or
%     string arrays and must have the same number of observations as y.
%
% Name-Value Options
%
%   cov: character vector specifying the method of estimating within-reader
%     rating covariances as 'DeLong', 'jackknife', or 'unbiased' [default:
%     'jackknife']
%
% Details
%
%   Readers and cases are treated as random factors by default. Either one
%   may be designated as fixed in calls to mrmc with the syntax
%   FixedVariate(<variable>), where <variable> is the reader or case id
%   variable.
%
% Return Value
%
%   MRMCFit class object of data that can be used to display and
%   statistically compare the reader performance metrics.
%
% Example
%
%   load VanDyke
%   y = ROCAUCVariate(VanDyke.truth, VanDyke.rating);
%   fit = mrmc(y, VanDyke.treatment, VanDyke.reader, VanDyke.case,...
%              'cov', 'unbiased');
%   plot(fit)
%   disp(fit)
%   summary(fit)
%
% See also: ROCAUCVariate, ROCEUVariate, ROCLRnegVariate, ROCLRposVariate,
% SensitivityVariate, SpecificityVariate, FixedVariate, MRMCFit

  arguments
    y PerformanceVariate
    test (:, 1)
    reader (:, 1)
    id (:, 1)
    options.cov char {...
      mustBeMember(options.cov, {'DeLong', 'jackknife', 'unbiased'})...
    } = 'jackknife'
  end

  fixed.reader = isa(reader, "FixedVariate");
  fixed.id = isa(id, "FixedVariate");
  assert(~(fixed.reader && fixed.id),...
         'Only one of reader or case may be fixed.')
  
  mrmc_factors = table(...
    categorical(test), categorical(reader), categorical(id),...
    'VariableNames', {'test', 'reader', 'id'}...
  );
  mrmc_factors = mrmc_factors(y.keep, :);
  design = get_design(mrmc_factors.test, mrmc_factors.reader, mrmc_factors.id);
  design.fixed = fixed;
  
  assert(~(design.id_in_reader && design.id_in_test), ...
         'Unsupported design with id nested within both reader and test')
  %assert(~design.reader_in_test,...
  %       'Unsupported design with reader nested within test.')

  data = table();
  [mrmc_factors.group, readerids, testids] = findgroups(...
    mrmc_factors.reader, mrmc_factors.test...
   );
  data.reader = readerids;
  data.test = testids;
  data.y = splitapply(y.metric, y.truth, y.rating, mrmc_factors.group);
  tbl = tabulate(mrmc_factors.group);
  data.N = tbl(:, 2);

  design.cov = options.cov;
  switch design.cov
    case 'DeLong'
      cov_func = @cov_DeLong;
    case 'jackknife'
      cov_func = @cov_jackknife;
    case 'unbiased'
      cov_func = @cov_unbiased;
  end
  V = cov_func(y, mrmc_factors.group, mrmc_factors.id, design.balanced3D);

  nested = [0 0; 0 0];
  if (design.reader_in_test)
    nested(1, :) = [0 1];
  end
  [~, tbl, anova.stats] = anovan(...
    data.y, {data.reader data.test},...
    'model', 'full', 'nested', nested, 'varnames', {'reader', 'test'},...
    'display', 'off'...
  );
  anova.tbl = cell2table(tbl(2:end, :), 'VariableNames', tbl(1, :));

  test_names = categories(data.test);
  n = length(test_names);
  testfits(1:n) = MRMCTestFit(char(), y, struct(), table(), struct(), nan);
  for i = 1:n
    test_name = test_names{i};
    inds = data.test == test_name;
    test_data = data(inds, {'y' 'reader'});
    [~, tbl, test_anova.stats] = anova1(test_data.y, test_data.reader, 'off');
    test_anova.tbl = cell2table(tbl(2:end, :), 'VariableNames', tbl(1, :));
    testfits(i) = MRMCTestFit(...
      test_name, y, design, test_data, test_anova, V(inds, inds)...
    );
  end

  fit = MRMCFit(y, mrmc_factors, design, data, anova, V, testfits);

end
