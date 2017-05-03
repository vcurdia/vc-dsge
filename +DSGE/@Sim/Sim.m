classdef Sim < handle
% DSGE.Sim class
% 
% See also:
% setupMyDSGE, DSGE.Model, DSGE.Prior, DSGE.Posterior
%
% ...........................................................................
%
% Created: May 2, 2017 by Vasco Curdia
% 
% Copyright (C) 2017 Vasco Curdia
    
    properties
        Model
        Prior
        Post
        Data
        Suffix
        Time2Show
        TimeLabels
    end
   
    methods
        function obj = Sim(model,prior,post,data)
            if nargin>0, obj.Model = model; end
            if nargin>1, obj.Prior = prior; end
            if nargin>2, obj.Post = post; end
            if nargin>3, obj.Data = data; end
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

