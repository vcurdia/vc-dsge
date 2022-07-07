classdef Param < matlab.mixin.Copyable

% DSGE.Param class
% 
% DSGE object representing parameters, used in the Model object.
%
% See also:
% DSGE.Model
%
% Created: November 7, 2016
% Copyright 2016-2018 Vasco Curdia
    
    properties 
% Names List of parameter names to be used in equations and commands.
        Names 

% PrettyNames - List of parameter names formatted for LaTeX
%  Use this property to create LaTeX output and figures.
        PrettyNames

% Values Vector of parameter values used in calibrated simulations.
        Values
        
% PriorDist String representing prior distribution type.
        PriorDist
        
% PriorMean Prior mean.
        PriorMean
        
% PriorSD Prior standard deviation.
        PriorSD
    end
    
    properties (SetAccess = protected)
        %N Number of parameters in instance
        N = 0;
        
        %NameLength maximum length of Names
        NameLength = 0;
    end
    
    methods
        
        function obj = Param(p)
            if nargin>0
                [np,nc] = size(p);
                obj.Names = p(:,1);
                if ismember(nc,[2,3])
                    obj.Values = [p{:,2}]';
                elseif nc>3
                    obj.PriorDist = p(:,2);
                    obj.PriorMean = [p{:,3}]';
                    for j=1:np
                        if isempty(p{j,4})
                            p{j,4} = 0;
                        end
                    end
                    obj.PriorSD = [p{:,4}]';
                    obj.Values = obj.PriorMean;
                end
                if ismember(nc,[3,5])
                    obj.PrettyNames = p(:,nc);
                end
            end
        end
        
        function obj = set.Names(obj,names)
            obj.Names = names;
            obj.N = length(names);
            obj.NameLength = max([cellfun('length',obj.Names)]);
            if length(obj.PrettyNames)~=obj.N
                obj.PrettyNames = names;
            end
        end
        
        function obj = set.PrettyNames(obj,prettynames)
            if length(prettynames)==obj.N
                obj.PrettyNames = prettynames;
            else
                error('Length of PrettyNames must match number of parameters.')
            end
        end
        
        function obj = set.Values(obj,values)
            if length(values)==obj.N
                obj.Values = values;
            else
                error('Length of Values must match number of parameters.')
            end
        end      
        
        function add(obj,p)
            if ~iscell(p)
                error('Parameters to add need to be in cell array.')
            end
            [np,nc] = size(p);
            names = [obj.Names;p(:,1)];
            if ismember(nc,[3,5])
                prettynames = [obj.PrettyNames;p(:,nc)];
            else
                prettynames = [obj.PrettyNames;p(:,1)];
            end
            values = obj.Values;
            obj.Names = names;
            obj.PrettyNames = prettynames;
            if ismember(nc,[2,3])
                obj.Values = [obj.Values;[p{:,2}]'];
            elseif nc>3
                obj.Values = [obj.Values;[p{:,3}]'];
                obj.PriorDist = [obj.PriorDist;p(:,2)];
                obj.PriorMean = [obj.PriorMean;[p{:,3}]'];
                for j=1:np
                    if isempty(p{j,4})
                        p{j,4} = 0;
                    end
                end
                obj.PriorSD = [obj.PriorSD;[p{:,4}]'];
            end
        end

        function showvalues(obj,x)
            if nargin<2
                x = obj.Values;
            end
            for jP=1:obj.N
                fprintf(['%',int2str(obj.NameLength),'s'],obj.Names{jP});
                fprintf('  %8.4f',x(jP,:));
                fprintf('\n');
            end
            fprintf('\n');
        end
        
        function x = get(obj,pname,x)
            if nargin<3, x = obj.Values; end
            x = x(ismember(obj.Names,pname),:);
        end

        function x = getvalues(obj,p)
        % loop allows for duplicates so that the last within the duplicates is
        % the one that remains.
        % invlaid param names ignored
            x = obj.Values;
            if nargin>1
                np = size(p,1);
                for jp=1:np
                    idx = ismember(obj.Names,p(jp,1));
                    if any(idx)
                        x(idx) = p{jp,2};
                    end
                end
            end
        end
    
        function p1 = merge(obj,varargin)
            p1 = copy(obj);
            isDuplicates = false;
            for j=1:nargin-1
                p = varargin{j};
                if ~strcmp(class(p),'DSGE.Param')
                    error('Cannot merge Param. Input needs to be instance of DSGE.Param')
                end
                if p.N==0, continue, end
                if p1.N==0
                    p1 = copy(p);
                    continue
                end
                tf = ismember(p.Names,p1.Names);
                isDuplicates = (isDuplicates || any(tf));
                if ~all(tf)
                    names = [p1.Names;p.Names(~tf)];
                    prettynames = [p1.PrettyNames;p.PrettyNames(~tf)];
                    values = [p1.Values;p.Values(~tf)];
                    priordist = [p1.PriorDist;p.PriorDist(~tf)];
                    priormean = [p1.PriorMean;p.PriorMean(~tf)];
                    priorsd = [p1.PriorSD;p.PriorSD(~tf)];
                    p1.Names = names;
                    p1.PrettyNames = prettynames;
                    p1.Values = values;
                    p1.PriorDist = priordist;
                    p1.PriorMean = priormean;
                    p1.PriorSD = priorsd;
                end
            end
%             if isDuplicates
%                 fprintf('Duplicate variable names found.\n')
%             end
        end
        
        function [p1,idx] = subset(obj,names)
            [tf,idx] = obj.ismember(names);
            idx = idx(tf);
            p1 = DSGE.Param(obj.Names(idx));
            p1.Values = obj.Values(idx);
            p1.PrettyNames = obj.PrettyNames(idx);
            p1.PriorDist = obj.PriorDist(idx);
            p1.PriorMean = obj.PriorMean(idx);
            p1.PriorSD = obj.PriorSD(idx);
        end
        
        function [tf,idx] = ismember(obj,names)
            if isa(names,'DSGE.Param'), names = names.Names; end
            [tf,idx] = ismember(names,obj.Names);
        end
        

    end %methods
    
end %class
