function FigPanels = GenVarFigPanels(obj,varargin)

% GenVarFigPanels
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
    op.FigShape = [1,1];
    op.PanelNames = {};
    op.PanelVar = {};
    
    %% Update options
    op = UpdateOptions(op,varargin{:});
    
    %% Check options
    if isempty(op.PanelNames)
        if obj.ObsVar.N>0, op.PanelNames{end+1} = 'ObsVar';end
        op.PanelNames{end+1} = 'StateVar';
        if obj.AuxVar.N>0, op.PanelNames{end+1} = 'AuxVar';end
    end
    nPanels = length(op.PanelNames);
    if isempty(op.PanelVar)
        for j=1:nPanels
            op.PanelVar{j} = obj.(op.PanelNames{j});
        end
    end
    
    %% Prepare panels
    FigPanels = struct;
    nVar = prod(op.FigShape);

    jP = 0;
    for jList=1:length(PanelList)
        Listj = PanelList{jList};
        nPanelj = ceil(obj.(Listj).N/obj.FigPanelMaxVar);
        for j=1:nPanelj
            jP = jP+1;
            if obj.FigPanelMaxVar==1
                obj.FigPanels(jP).Var = obj.(Listj).Names(j);
                obj.FigPanels(jP).PrettyNames = obj.(Listj).PrettyNames(j);
                obj.FigPanels(jP).Title = sprintf('%s_%s',Listj,...
                                                  obj.FigPanels(jP).Var{1});
            else
                if nPanelj>1
                    obj.FigPanels(jP).Title = sprintf('%s_%.0f',Listj,j);
                else
                    obj.FigPanels(jP).Title = sprintf('%s',Listj);
                end
                obj.FigPanels(jP).Var = obj.(Listj).Names(...
                    (j-1)*obj.FigPanelMaxVar+1:min(j*obj.FigPanelMaxVar,...
                                                   obj.(Listj).N));
                obj.FigPanels(jP).PrettyNames = obj.(Listj).PrettyNames(...
                    (j-1)*obj.FigPanelMaxVar+1:min(j*obj.FigPanelMaxVar,...
                                                   obj.(Listj).N));
            end
        end
    end





end

