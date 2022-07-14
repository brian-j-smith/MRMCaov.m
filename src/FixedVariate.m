classdef FixedVariate
  
  properties
    value categorical
  end
  
  methods
    
    function obj = FixedVariate(x)

% Description
%
% Designate readers or cases as fixed factors in an MRMC analysis.
%
% Usage
%
%   obj = FixedVariate(x)
%
% Input Arguments
%
%   x: column array of reader or case true binary statuses.
%
%   rating: column array of the grouping variable for readers or cases.
%
% Return Value
%
%   FixedVariate class object.
%
% See also: mrmc

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
