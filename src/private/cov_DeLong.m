function V = cov_DeLong(y, group, ~, balanced)

  arguments
    y ROCAUCVariate
    group (:, 1) double
    ~
    balanced (1, 1) logical
  end
  
  assert(~length(y.partial), 'DeLong covariance undefined for partial AUC.')
  assert(balanced, 'DeLong covariance requires a balanced study design.')

  comps = splitapply(@varcomp_Sen, y.truth, y.rating, group);
  
  V10 = horzcat(comps.v10);
  n_pos = size(V10, 1);
  V10 = (V10' * V10) / (n_pos - 1);
  
  V01 = horzcat(comps.v01);
  n_neg = size(V01, 1);
  V01 = (V01' * V01) / (n_neg - 1);
  
  V = V10 / n_pos + V01 / n_neg;
  
end


function comps = varcomp_Sen(truth, rating)

  is_pos = truth == category(truth, 2);
  rating_pos = rating(is_pos);
  rating_neg = rating(~is_pos);
  
  inds = combvec(1:length(rating_pos), 1:length(rating_neg))';
  inds_pos = inds(:, 1);
  inds_neg = inds(:, 2);
  scores = c_score(rating_pos(inds_pos), rating_neg(inds_neg));
  
  auc = roc_auc(truth, rating);
  
  v10 = splitapply(@mean, scores, inds_pos) - auc;
  v01 = splitapply(@mean, scores, inds_neg) - auc;
  
  comps.v10 = v10;
  comps.v01 = v01;
  
end
