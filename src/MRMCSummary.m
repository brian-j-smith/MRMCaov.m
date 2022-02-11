classdef MRMCSummary
  
  properties
    response char
    design struct
    alpha double
    vcov table
    test_equality table
    test_means table
    test_diffs table
    reader_test_diffs table
    reader_means table
  end
  
  methods
    
    function obj = MRMCSummary(fit, alpha, vcov, test_equality,...
                               test_means, test_diffs,...
                               reader_test_diffs, reader_means)
      obj.response = class(fit.y);
      obj.design = fit.design;
      obj.alpha = alpha;
      obj.vcov = vcov;
      obj.test_equality = test_equality;
      obj.test_means = test_means;
      obj.test_diffs = test_diffs;
      obj.reader_test_diffs = reader_test_diffs;
      obj.reader_means = reader_means;
    end
    
    
    function disp(x)

      fprintf('Multi-Reader Multi-Case Analysis of Variance\n')
      if x.design.id_in_reader
        str = 'cases nested within readers';
      elseif x.design.id_in_test
        str = 'cases nested within tests';
      else
        str = 'factorial';
      end
      fprintf('Experimental design: %s\n', str)
      types = {'random' 'fixed'};
      fprintf('Factor types: %s readers and %s cases\n',...
              types{x.design.fixed.reader + 1},...
              types{x.design.fixed.id + 1})
      fprintf('Response: %s\n', x.response)
      fprintf('Covariance method: %s\n', x.design.cov)
      fprintf('Confidence interval level: %s%%\n\n',...
              string(100 * (1 - x.alpha)))
      fprintf(['Obuchowski-Rockette variance component and covariance '...
               'estimates:\n\n'])
      if x.design.fixed.id
        fprintf('Not applicable because cases are fixed\n\n')
      else
        disp(x.vcov)
      end
      
      fprintf('ANOVA global statistical test of equal tests:\n\n')
      disp(x.test_equality)

      fprintf('Pairwise test differences:\n\n')
      disp(x.test_diffs)

      fprintf('Test means based only on the data for each one:\n\n')
      disp(x.test_means)

      if ~isempty(x.reader_test_diffs)
        fprintf('Reader-specific pairwise test differences:\n\n')
        disp(x.reader_test_diffs)
      end

      if ~isempty(x.reader_means)
        fprintf('Reader-specific test results:\n\n')
        disp(x.reader_means)
      end

    end
    
  end
  
end