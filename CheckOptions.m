function dsge = CheckOptions(dsge)

% SetDefaultOptions
%
% Checks options and if needed initializes them.
%
% Note: Every time that some option is changed that can influence multiple 
%       actions (e.g. Panels used in simulations) it is advisable to run this 
%       function again.
%
% See also:
% vcDSGE
%
% ...........................................................................
%
% Created: February 24, 2016 by Vasco Curdia
% 
% Copyright (C) 2016 Vasco Curdia

%% -------------------------------------------------------------------

%% Options for Simulations

if isfield(dsge.Options,'Sim')
    op = dsge.Options.Sim;
else
    op = struct; 
end 

if ~isfield(op,'UseDist'), op.UseDist = 'PriorPrc500'; end
if ~isfield(op,'nDraws')
    if ismember(op.UseDist,{'PriorDraws','PostDraws'})
        op.nDraws = 1000;
    else
        op.nDraws = 1;
    end
end

if ~isfield(op,'Bands2Show'), op.Bands2Show = [50,70,90]; end

if ~isfield(op,'PanelList')
    op.PanelList = {'ObsVar','StateVar','AuxVar'};
end
if ~isfield(op,'PanelMaxVar'), op.PanelMaxVar = 9; end
%if ~isfield(op,'PanelFigShape'), op.PanelFigShape = [3,3]; end
if ~isfield(op,'Panels')
    jP = 0;
    for jList=1:length(op.PanelList)
        Listj = op.PanelList{jList};
        for j=(1:ceil(dsge.n.(Listj)/op.PanelMaxVar))
            jP = jP+1;
            if op.PanelMaxVar==1
                op.Panels(jP).Var = dsge.(Listj).Names(j);
                op.Panels(jP).PrettyNames = dsge.(Listj).PrettyNames(j);
                op.Panels(jP).Title = sprintf('%s_%s',Listj,...
                                              op.Panels(jP).Var{1});
%                op.Panels(jP).FigShape = [1,1];
            else
                op.Panels(jP).Title = sprintf('%s_%.0f',Listj,j);
                op.Panels(jP).Var = dsge.(Listj).Names(...
                    (j-1)*op.PanelMaxVar+1:min(j*op.PanelMaxVar,...
                                               dsge.n.(Listj)));
                op.Panels(jP).PrettyNames = dsge.(Listj).PrettyNames(...
                    (j-1)*op.PanelMaxVar+1:min(j*op.PanelMaxVar,...
                                               dsge.n.(Listj)));
%                op.Panels(jP).FigShape = op.PanelFigShape;
            end
        end
    end
end
op.nPanels = length(op.Panels);
for jP=1:op.nPanels
    op.Panels(jP).nVar = length(op.Panels(jP).Var);
    if ~isfield(op.Panels(jP),'PrettyNames')
        op.Panels(jP).PrettyNames = op.Panels(jP).Var;
    end
    if ~isfield(op.Panels(jP),'Scale') || isempty(op.Panels(jP).Scale)
        op.Panels(jP).Scale = ones(1,op.Panels(jP).nVar);
    end
end

if ~isfield(op,'Shocks2Show'), op.Shocks2Show = dsge.ShockVar.Names; end
op.nShocks2Show = length(op.Shocks2Show);

if ~isfield(op,'VarScale'), op.VarScale = 1; end

dsge.Options.Sim = op;


%% Figures

if isfield(dsge.Options,'Fig')
    op = dsge.Options.Fig;
else
    op = struct; 
end 

if ~isfield(op,'Visible'), op.Visible = 'off'; end
if ~isfield(op,'YSlack'), op.YSlack = 0.05; end
if ~isfield(op,'YMinScale'), op.YMinScale = 0; end
if ~isfield(op,'KeepEPS'), op.KeepEPS = 0; end
if ~isfield(op,'OpenPDF'), op.OpenPDF = 0; end

dsge.Options.Fig = op;


%% Parameter Tables

if isfield(dsge.Options,'ParTable')
    op = dsge.Options.ParTable;
else
    op = struct; 
end 

if ~isfield(op,'MoveLeft'), op.MoveLeft = 1; end
if ~isfield(op,'Precision'), op.Precision = 3; end
if ~isfield(op,'MaxRows'), op.MaxRows = 35; end
if ~isfield(op,'Lines'), op.Lines = []; end

dsge.Options.ParTable = op;

%% -------------------------------------------------------------------
