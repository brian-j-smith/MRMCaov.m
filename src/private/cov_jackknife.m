function V = cov_jackknife(y, group, id, ~)

  arguments
    y PerformanceVariate
    group (:, 1) double
    id (:, 1) categorical
    ~
  end

  id_names = categories(id);
  m = length(id_names);
  n = length(unique(group));
  X = NaN(m, n);
  wb = waitbar(0, 'Computing jackknife covariance matrix');
  for i = 1:m
    keep = id ~= id_names{i};
    X(i, :) = splitapply(y.metric, y.truth(keep), y.rating(keep), group(keep));
    waitbar(i / m, wb);
  end
  close(wb);
  X = X - mean(X);
  V = ((m - 1) / m) * (X' * X);
  
end
