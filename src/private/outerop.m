function X = outerop(x, y, func)

  arguments
    x (:, 1)
    y (:, 1)
    func function_handle
  end
  
  X = func(repmat(x, 1, length(y)), y');
  
end
