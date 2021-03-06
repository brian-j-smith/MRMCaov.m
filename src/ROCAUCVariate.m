classdef ROCAUCVariate < PerformanceVariate
  
  properties (SetAccess = private)
    partial char {...
      mustBeMember(partial, {'', 'sensitivity', 'specificity'})...
    } = ''
    min (1,1) double {mustBeInRange(min, 0, 1)}
    max (1,1) double {mustBeInRange(max, 0, 1)}
    normalize (1,1) logical
  end
  
  methods
    
    function obj = ROCAUCVariate(truth, rating, options)

      arguments
        truth
        rating
        options.partial = ''
        options.min = 0
        options.max = 1
        options.normalize = false
      end

% Description
%
% Designates area under the ROC curve as the reader performance metric to
% be computed from true case statuses and reader ratings for an MRMC
% analysis.
%
% Usage
%
%   obj = ROCAUCVariate(truth, rating, options)
%
% Input Arguments
%
%   truth: column array of true binary statuses.
%
%   rating: column array of numeric ratings.
%
% Name-Value Options
%
%   partial: character vector specifying whether to compute area under the
%     entire curve (''), over a range of sensitivities ('sensitivity'), or
%     over a ranges of specificities ('specificity') [default: ''].
%
%   min: numeric value from 0 to 1 for the minimum sensitivity or
%     specificity if computing partial area [default: 0].
%
%   max: numeric value from 0 to 1 for the maximum sensitivity or
%     specificity if computing partial area [default: 1].
%
%   normalize: logical indicating whether to divide partial area by
%     max - min [default: false].
%
% Return Value
%
%   ROCAUCVariate class object.
%
% See also: mrmc
      metric = @(truth, rating) roc_auc(...
        truth, rating, 'partial', options.partial, 'min', options.min,...
        'max', options.max, 'normalize', options.normalize);
      obj@PerformanceVariate(categorical(truth), double(rating), metric)
      assert(length(categories(obj.truth)) == 2,...
             'Argument ''truth'' must have 2 categories.')
      obj.label = 'ROC AUC';
      obj.partial = options.partial;
      obj.min = options.min;
      obj.max = options.max;
      obj.normalize = options.normalize;
    end
    
    
    function plot(obj, groups)
      
      arguments
        obj
        groups table = table()
      end
      
      if width(groups)
        group1 = categorical(groups{:, 1});
      else
        group1 = categorical(ones(length(obj), 1));
      end
      
      if width(groups) > 1
        [group2, group2ids] = findgroups(groups(:, 2:end));
        group2 = categorical(group2);
        totplots = length(categories(group2));
        nplots = ceil(sqrt(totplots));
        mplots = ceil(totplots / nplots);
      else
        group2 = categorical(ones(length(obj), 1));
        group2ids = '1';
        totplots = 1;
        mplots = 1;
        nplots = 1;
      end
      
      for i = 1:totplots
        
        inds2 = group2 == category(group2, i);
        
        subplot(mplots, nplots, i)

        inds1 = group1 == category(group1, 1) & inds2;
        truth = obj.truth(inds1);
        rating = obj.rating(inds1);
        [fpr, tpr] = perfcurve2(truth, rating, category(truth, 2));
        
        plot(fpr, tpr)
        subplot_title(group2ids, i)
        xlabel({'False Positive Rate'})
        ylabel({'True Positive Rate'})
        axis equal
        pad = 0.05 * [-1 1 -1 1];
        axis([0 1 0 1] + pad)

        hold on
        for j = 2:length(categories(group1))
          inds1 = group1 == category(group1, j) & inds2;
          truth = obj.truth(inds1);
          rating = obj.rating(inds1);
          [fpr, tpr] = perfcurve2(truth, rating, category(truth, 2));
          
          subplot_title(group2ids, i)
          plot(fpr, tpr)
        end
        hold off

        if width(groups)
          lgd = legend(categories(group1), 'Location', 'southeast');
          title(lgd, groups.Properties.VariableNames{1})
        end
    
      end
      
    end
    
    
    function value = trunc(~, x)
      value = max(0, min(x, 1));
    end
    
  end
  
end


function [X, Y] = perfcurve2(labels, scores,posclass)
  if length(unique(labels)) == 2
    [X, Y] = perfcurve(labels, scores, posclass);
  else
    X = [0.5];
    Y = [0.5];
  end
end


function subplot_title(groupids, i)
  if istable(groupids)
    names = groupids.Properties.VariableNames;
    values = groupids{i, :};
    title(strjoin(strcat(names, " = ", string(values)), {', '}))
  end
end
