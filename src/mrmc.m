function fit = mrmc(y, test, reader, id, options)

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
  is_fixed = fixed.reader || fixed.id;
  assert(~(fixed.reader && fixed.id),...
         'Only one of reader or case may be fixed.')
  
  mrmc_factors = table(...
    categorical(test), categorical(reader), categorical(id),...
    'VariableNames', {'test', 'reader', 'id'}...
  );
  mrmc_factors = mrmc_factors(y.keep, :);
  design = get_design(mrmc_factors.test, mrmc_factors.reader, mrmc_factors.id);
  design.fixed = fixed;
  
  is_nested = design.id_in_reader || design.id_in_test || design.reader_in_test;
  assert(~is_nested, ['Support for nested designs is under development and '...
                      'not currently available.'])

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

  [~, tbl, anova.stats] = anovan(...
    data.y, {data.reader data.test},...
    'model', 'interaction', 'varnames', {'reader', 'test'}, 'display', 'off'...
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
