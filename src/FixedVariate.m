classdef FixedVariate
  
  properties
    value categorical
  end
  
  methods
    
    function obj = FixedVariate(x)
      obj.value = categorical(x);
    end
    
    function disp(x)
      disp(x.value)
    end
    
    function B = categorical(A)
      B = A.value;
    end
    
    function sz = size(x, varargin)
      sz = size(x.value, varargin{:});
    end
    
  end
  
end
