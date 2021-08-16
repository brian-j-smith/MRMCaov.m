function value = roc_auc(truth, rating, options)

  arguments
    truth (:, 1)
    rating(:, 1) double
    options.partial char {...
      mustBeMember(options.partial, {'', 'sensitivity', 'specificity'})...
    } = ''
    options.min (1,1) double {mustBeInRange(options.min, 0, 1)} = 0
    options.max (1,1) double {mustBeInRange(options.max, 0, 1)} = 1
    options.normalize (1,1) logical = false
  end
  
  truth = categorical(truth);
  partial = options.partial;
  min = options.min;
  max = options.max;
  normalize = options.normalize;
  
  assert(max > min, 'Value of ''max'' must be greater than ''min''.')
  
  [fpr, tpr] = perfcurve(truth, rating, category(truth, 2));
  
  switch partial
    case 'sensitivity'
      [y, x] = rocsubset(tpr, fpr, min, max);
    case 'specificity'
      [x, y] = rocsubset(fpr, tpr, 1 - max, 1 - min);
    otherwise
      min = 0;
      max = 1;
      x = fpr;
      y = tpr;
  end
  
  value = trapz(x, y);
  if normalize
    value = value / (max - min);
  end
  
end


function [xsub, ysub] = rocsubset(x, y, xmin, xmax)
  keep = x >= xmin & x <= xmax;
  xsub = x(keep);
  ysub = y(keep);
  if (~any(x == xmin))
    xsub = [xmin; xsub];
    ysub = [rocinterp(x, y, xmin); ysub];
  end
  if (~any(x == xmax))
    xsub = [xsub; xmax];
    ysub = [ysub; rocinterp(x, y, xmax)];
  end
end


function ynew = rocinterp(x, y, xnew)

  inds = find(x <= xnew);
  ind1 = inds(end);
  x1 = x(ind1);
  
  inds = find(x >= xnew);
  ind2 = inds(1);
  x2 = x(ind2);

  ynew = ((x2 - xnew) * y(ind1) + (xnew - x1) * y(ind2)) / (x2 - x1);
  
end
