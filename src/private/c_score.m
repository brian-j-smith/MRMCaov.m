function value = c_score(x_pos, x_neg)
  value = (x_pos > x_neg) + 0.5 * (x_pos == x_neg);
end
