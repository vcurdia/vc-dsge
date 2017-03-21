function FigPanels = setvarfigpanels(obj,varargin)

% setvarfigpanels
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
op.FigShape = {};
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
if isempty(op.FigShape), op.FigShape = cell(1,nList); end
if ~iscell(op.FigShape)
    FigShape = op.FigShape;
    op.FigShape = cell(1,nList);
    for j=1:nList, op.FigShape{j} = FigShape; end
    clear FigShape
end
if ~iscell(op.Scale)
    Scale = op.Scale;
    op.Scale = cell(1,nList);
    for j=1:nList, op.Scale{j} = Scale; end
    clear Scale
end
for j=1:nList
    if ~isempty(op.FigShape{j}), continue, end
    nVarj = obj.(op.PanelList{j}).N;
    if nVarj==1
        op.FigShape{j} = [1,1];
    elseif nVarj<=4*2
        op.FigShape{j} = [2,2];
    elseif nVarj<=9*3
        op.FigShape{j} = [3,3];
    else
        op.FigShape{j} = [4,4];
    end
end

%% Prepare panels
FigPanels = struct;
jP = 0;
for jL=1:nList
    Lj = op.PanelList{jL};
    nMaxVar = prod(op.FigShape{jL});
    nPj = ceil(obj.(Lj).N/nMaxVar);
    for j=1:nPj
        jP = jP+1;
        if nMaxVar==1
            FigPanels(jP).Title = sprintf('%s_%s',Lj,obj.(Lj).Names{j});
            FigPanels(jP).Names = obj.(Lj).Names(j);
            FigPanels(jP).PrettyNames = obj.(Lj).PrettyNames(j);
        else
            if nPj>1
                FigPanels(jP).Title = sprintf('%s_%.0f',Lj,j);
            else
                FigPanels(jP).Title = sprintf('%s',Lj);
            end
            FigPanels(jP).Names = obj.(Lj).Names(...
                (j-1)*nMaxVar+1:min(j*nMaxVar,obj.(Lj).N));
            FigPanels(jP).PrettyNames = obj.(Lj).PrettyNames(...
                (j-1)*nMaxVar+1:min(j*nMaxVar,obj.(Lj).N));
        end
        FigPanels(jP).FigShape = op.FigShape{jL};
        FigPanels(jP).Scale = ...
            repmat(op.Scale{jL},1,length(FigPanels(jP).Names));
    end
end

