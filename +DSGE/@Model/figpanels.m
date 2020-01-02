function FigPanels = figpanels(obj,varargin)

% figpanels
% 
% Generates figure panels to show variables from DSGE model
% 
% See also:
% 
%
% ...........................................................................
% 
% Created: January 12, 2017 by Vasco Curdia
% 
% Copyright 2017 by Vasco Curdia

%% Options
op.PanelList = {};
op.Shape = {3,2};
op.Scale = 1;

%% Update options
op = updateoptions(op,varargin{:});

%% Check options
if isempty(op.PanelList)
    if obj.ObsVar.N>0, op.PanelList{end+1} = 'ObsVar'; end
    op.PanelList{end+1} = 'StateVar';
    if obj.AuxVar.N>0, op.PanelList{end+1} = 'AuxVar';end
end
nList = length(op.PanelList);
if isempty(op.Shape), op.Shape = cell(1,nList); end
if ~isempty(op.Shape{1}) && ~iscell(op.Shape{1})
    Shape = op.Shape;
    op.Shape = cell(1,nList);
    for j=1:nList, op.Shape{j} = Shape; end
    clear Shape
end
if ~iscell(op.Scale)
    Scale = op.Scale;
    op.Scale = cell(1,nList);
    for j=1:nList, op.Scale{j} = Scale; end
    clear Scale
end
for j=1:nList
    if ~isempty(op.Shape{j}), continue, end
    nVarj = obj.(op.PanelList{j}).N;
    if nVarj==1
        op.Shape{j} = {1,1};
    elseif nVarj<=4*2
        op.Shape{j} = {2,2};
    elseif nVarj<=9*3
        op.Shape{j} = {3,3};
    else
        op.Shape{j} = {4,4};
    end
end

%% Prepare panels
FigPanels = struct;
jP = 0;
for jL=1:nList
    Lj = op.PanelList{jL};
    nMaxVar = prod([op.Shape{jL}{:}]);
    nPj = ceil(obj.(Lj).N/nMaxVar);
    for j=1:nPj
        jP = jP+1;
        if nMaxVar==1
            FigPanels(jP).Title = sprintf('%s-%s',lower(Lj),obj.(Lj).Names{j});
            FigPanels(jP).Names = obj.(Lj).Names(j);
            FigPanels(jP).PrettyNames = obj.(Lj).PrettyNames(j);
        else
            if nPj>1
                FigPanels(jP).Title = sprintf('%s-%.0f',lower(Lj),j);
            else
                FigPanels(jP).Title = sprintf('%s',lower(Lj));
            end
            FigPanels(jP).Names = obj.(Lj).Names(...
                (j-1)*nMaxVar+1:min(j*nMaxVar,obj.(Lj).N));
            FigPanels(jP).PrettyNames = obj.(Lj).PrettyNames(...
                (j-1)*nMaxVar+1:min(j*nMaxVar,obj.(Lj).N));
        end
        FigPanels(jP).Shape = op.Shape{jL};
        FigPanels(jP).Scale = ...
            repmat(op.Scale{jL},1,length(FigPanels(jP).Names));
    end
end

