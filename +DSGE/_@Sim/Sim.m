classdef Sim < handle

% DSGE.Sim class
% 
% Unfinished class to implement model simulations as an object.
%
% See also:
% setupMyDSGE, DSGE.Model, DSGE.Prior, DSGE.Posterior, DSGE.Data
%
% Created: May 2, 2017
% Copyright 2017 Vasco Curdia
    
    properties
        Name
        Model
        Prior
        Post
        Data
        Dist
        NDraws = 1;
        Time2Show
        TimeLabels
        IRF = 1;
        VD = 1;
        States = 1;
        SD = 1;
        IRFNSteps = 25;
        IRFTickStep = 4;
        IRFFigPanels
        Shocks2Show
        ShockSize
        IRFPlotDir = 'Plots_IRF';
        Fig = figdefaultoptions;
    end
   
    properties (SetAccess = protected)
        TimeTracker
    end
    
    methods
        function obj = Sim(model,prior,post,data)
            obj.TimeTracker = TimeTracker;
            if nargin>0 && ~isempty(model), obj.Model = model; end
            if nargin>1 && ~isempty(prior), obj.Prior = prior; end
            if nargin>2 && ~isempty(post), obj.Post = post; end
            if nargin>3 && ~isempty(data), obj.Data = data; end
        end
        
        function set.Model(obj,m)
            obj.Model = m;
            if isempty(obj.Name), obj.Name = m.Name; end
        end
        
        function set.Prior(obj,p)
            obj.Prior = p;
            obj.checkparam('Prior')
        end
        
        function set.Post(obj,p)
            obj.Post = p;
            obj.checkparam('Post')
        end
        
        function set.Data(obj,data)
            obj.Data = data;
            obj.Data.Var = obj.Model.ObsVar.Names;
            if isempty(obj.Time2Show)
                obj.Time2Show = data.TimeIdx([1,end]);
            end
        end
        
        function checkparam(obj,pp)
            if isempty(obj.Model)
                error(['Cannot check parameter consistency without specifying ' ...
                       'model.'])
            end
            pm = obj.Model.Param; 
            p = obj.(pp).Model.Param;
            tf = (p.N==pm.N);
            names = {pm.Names{:},p.Names{:}};
            nl = max([cellfun('length',names)]);
            if tf, tf = strcmp(pm.Names,p.Names); end
            if ~all(tf==1)
                fprintf(['%4s',repmat([' %',int2str(nl),'s'],1,2)],...
                        '','model',pp)
                 for j=1:max(pm.N,p.N)
                    if j<=pm.N,pmj=pm.Names{j};else pmj = '';end
                    if j<=p.N,pj=p.Names{j};else pj = '';end
                    fprintf(['%3.0f:',repmat([' %',int2str(nl),'s'],1,2)],...
                            '',pmj,pj)
                end
                error('Parameters do not match. Cannot proceed.')
            end
        end
        
        function list(obj,slist)
            simList = {'IRF','VD','States','SD'};
            tf = ismember(simList,sList);
            for j=1:length(simList)
                obj.(simList{j}) = tf(j);
            end
        end
        
        function Fig = figdefaultoptions
            Fig = struct;
            Fig.Visible = 'off';
            Fig.Plot.LineWidth = 1.5;
        end
        
        function new = copy(obj)
            new = DSGE.Sim;
            % Copy all non-hidden properties.
            pList = properties(obj);
            for j = 1:length(pList)
                new.(pList{j}) = obj.(pList{j});
            end
        end
    
    end %methods
    
end %class

