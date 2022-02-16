function V = cov_unbiased(y, group, id, balanced)

  arguments
    y ROCAUCVariate
    group (:, 1) double
    id (:, 1) categorical
    balanced (1, 1) logical
  end
  
  assert(~length(y.partial), 'Unbiased covariance undefined for partial AUC.')
  
  if balanced
    func = @cov_unbiased_balanced;
  else
    func = @cov_unbiased_default;
  end
  V = func(y.truth, y.rating, group, id);
  
end


function V = cov_unbiased_default(truth, rating, group, id)

  n = length(unique(group));
  is_pos = truth == category(truth, 2);

  V = zeros(n);
  x = cell(n, 1);
  for i = 1:n
    x{i} = get_scores(rating, group == i, is_pos, id);
    for j = 1:i
      V(i, j) = get_cov(x{i}, x{j});
      V(j, i) = V(i, j);
    end
  end

end


function res = get_scores(rating, is_group, is_pos, id)
  is_group_pos = is_group & is_pos;
  is_group_neg = is_group & ~is_pos;
  res.id_pos = id(is_group_pos);
  res.id_neg = id(is_group_neg);
  res.scores = outerop(rating(is_group_pos), rating(is_group_neg), @c_score);
  res.sum_scores_pos = sum(res.scores, 2);
  res.sum_scores_neg = sum(res.scores, 1);
  res.sum_scores = sum(res.sum_scores_neg);
end


function res = get_cov(x, y)

  id_pos = intersect(x.id_pos, y.id_pos);
  id_neg = intersect(x.id_neg, y.id_neg);
  nxny = numel(x.scores) * numel(y.scores);
  n_pos = length(id_pos);
  n_neg = length(id_neg);

  delta = nxny / (...
    nxny -...
    size(x.scores, 1) * n_neg * size(y.scores, 1) -...
    size(x.scores, 2) * n_pos * size(y.scores, 2) +...
    n_pos * n_neg...
  );

  [~, inds_pos] = ismember(id_pos, x.id_pos);
  [~, inds_neg] = ismember(id_neg, x.id_neg);
  x.sum_scores_pos = x.sum_scores_pos(inds_pos);
  x.sum_scores_neg = x.sum_scores_neg(inds_neg);
  x.scores = x.scores(inds_pos, inds_neg);
  
  [~, inds_pos] = ismember(id_pos, y.id_pos);
  [~, inds_neg] = ismember(id_neg, y.id_neg);
  y.sum_scores_pos = y.sum_scores_pos(inds_pos);
  y.sum_scores_neg = y.sum_scores_neg(inds_neg);
  y.scores = y.scores(inds_pos, inds_neg);
  
  res = (...
    (1 - delta) * (x.sum_scores * y.sum_scores) +...
    delta * (sum(x.sum_scores_pos .* y.sum_scores_pos) +...
             sum(x.sum_scores_neg .* y.sum_scores_neg) -...
             sum(x.scores .* y.scores, "all"))...
  ) / nxny;

end


function V = cov_unbiased_balanced(truth, rating, group, id)

  [~, sort_inds] = sortrows(table(group, id));
  
  m = length(categories(id));
  n = length(rating) / m;
  rating_mat = reshape(rating(sort_inds), m, n);
  truth = truth(sort_inds(1:m));
  
  is_pos = truth == category(truth, 2);
  rating_mat_pos = rating_mat(is_pos, :);
  rating_mat_neg = rating_mat(~is_pos, :);
  n_pos = size(rating_mat_pos, 1);
  n_neg = size(rating_mat_neg, 1);
  
  scores = zeros(n_pos * n_neg, n);
  sum_scores_pos = zeros(n_pos, n);
  sum_scores_neg = zeros(n_neg, n);
  
  for i = 1:n
    x = outerop(rating_mat_pos(:, i), rating_mat_neg(:, i), @c_score);
    scores(:, i) = x(:);
    sum_scores_pos(:, i) = sum(x, 2);
    sum_scores_neg(:, i) = sum(x, 1);
  end
  
  delta = (n_neg * n_pos)^2 / (n_neg * (n_neg - 1) * n_pos * (n_pos - 1));
  sum_scores = sum(scores, 1);
  V = (...
    (1 - delta) * (sum_scores' * sum_scores) +...
    delta * (sum_scores_pos' * sum_scores_pos +...
             sum_scores_neg' * sum_scores_neg -...
             scores' * scores)...
  ) / (n_neg * n_pos)^2;
  
end
