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
    func = @cov_unbiased_unbalanced;
  end
  V = func(y.truth, y.rating, group, id);
  
end


function V = cov_unbiased_balanced(truth, rating, group, id)

  [~, sort_inds] = sortrows(table(group, id));
  
  m = length(categories(id));
  n = length(rating) / m;
  X = reshape(rating(sort_inds), m, n);
  truth = truth(sort_inds(1:m));
  
  is_pos = truth == category(truth, 2);
  X_pos = X(is_pos, :);
  X_neg = X(~is_pos, :);
  n_pos = size(X_pos, 1);
  n_neg = size(X_neg, 1);
  
  V = NaN(n);
  scores = cell(n, 1);
  for i = 1:n
    x = outerop(X_pos(:, i), X_neg(:, i), @c_score);
    scores{i} = x;
    V(i, i) = var_unbiased(x, n_pos, n_neg);
  end
  
  for i = 1:n
    for j = 1:(i - 1)
      x = scores{i} - scores{j};
      var_diff = var_unbiased(x, n_pos, n_neg);
      V(i, j) = (V(i, i) + V(j, j) - var_diff) / 2;
      V(j, i) = V(i, j);
    end
  end
  
end


function v = var_unbiased(x, n_pos, n_neg)
    q1 = sum(x.^2, "all");
    m1 = q1 / (n_pos * n_neg);
    q2 = sum(sum(x, 1).^2);
    m2 = (q2 - q1) / (n_pos * n_neg * (n_neg - 1));
    q3 = sum(sum(x, 2).^2);
    m3 = (q3 - q1) / (n_neg * n_pos * (n_pos - 1));
    q4 = sum(x, "all").^2;
    m4 = (q4 - q2 - q3 + q1) / (n_neg * n_pos * (n_neg - 1) * (n_pos - 1));
    v = m1 / (n_neg * n_pos) + m2 * (n_neg - 1) / (n_neg * n_pos) +...
      m3 * (n_pos - 1) / (n_neg * n_pos) +...
      m4 * ((n_neg - 1) * (n_pos - 1) / (n_neg * n_pos) - 1);
end


function V = cov_unbiased_unbalanced(truth, rating, group, id)

  data = table(rating, group, id,...
               'VariableNames', {'rating', 'group', 'id'});
  n = length(unique(group));
  
  is_pos = truth == category(truth, 2);
  data_pos = data(is_pos, :);
  data_neg = data(~is_pos, :);
  n_pos = length(unique(data_pos.id));
  n_neg = length(unique(data_neg.id));
  
  scores = cell(n, 1);
  for i = 1:n
    scores{i} = group_c_scores(data_pos, data_neg, i);
  end
  
  A1 = NaN(n);
  A2 = NaN(n);
  A3 = NaN(n);
  A4 = NaN(n);

  iter = 0;
  niter = n * (n + 1) / 2;
  wb = waitbar(0, 'Computing unbiased covariance matrix');
  for i = 1:n
    for j = 1:i
      X = scores{i}.value * scores{j}.value';
      same_id_pos = outerop(scores{i}.id_pos, scores{j}.id_pos, @eq);
      same_id_neg = outerop(scores{i}.id_neg, scores{j}.id_neg, @eq);
      
      A1(i, j) = mean(X(same_id_pos & same_id_neg));
      A1(j, i) = A1(i, j);
      A2(i, j) = mean(X(same_id_pos & ~same_id_neg));
      A2(j, i) = A2(i, j);
      A3(i, j) = mean(X(~same_id_pos & same_id_neg));
      A3(j, i) = A3(i, j);
      A4(i, j) = mean(X(~same_id_pos & ~same_id_neg));
      A4(j, i) = A4(i, j);
      
      iter = iter + 1;
      waitbar(iter / niter, wb);
    end
  end
  close(wb);
  
  V = (A1 + (n_neg - 1) * A2 + (n_pos - 1) * A3 +...
       (1 - n_pos - n_neg) * A4) / (n_pos * n_neg);
  
end


function scores = group_c_scores(data_pos, data_neg, which)

  is_group_pos = data_pos.group == which;
  rating_pos = data_pos.rating(is_group_pos);
  id_pos = data_pos.id(is_group_pos);
  
  is_group_neg = data_neg.group == which;
  rating_neg = data_neg.rating(is_group_neg);
  id_neg = data_neg.id(is_group_neg);
  
  inds = combvec(1:length(rating_pos), 1:length(rating_neg))';
  inds_pos = inds(:, 1);
  inds_neg = inds(:, 2);
  
  scores.value = c_score(rating_pos(inds_pos), rating_neg(inds_neg));
  scores.id_pos = id_pos(inds_pos);
  scores.id_neg = id_neg(inds_neg);
  
end
