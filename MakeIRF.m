function s = MakeIRF(s)

% MakeIRF
% 
% Generates IRF for DSGE model
% 
% See also:
% MakeIRFFcn
%
% .............................................................................
% 
% Created: February 9, 2016 by Vasco Curdia
% 
% Copyright 2016 by Vasco Curdia

%-----------------------------------------------------------------------------%

%% Preamble

Action = 'IRF';

fprintf('\n*** Making IRF\n')

% Set Timer
s.TimeElapsed.(Action) = toc();

%% Options

if isfield(s.Options,'IRF')
    op = s.Options.IRF;
else
    op = struct; 
end

if ~isfield(op,'UseDist'), op.UseDist = PriorPrc500; end
if ~isfield(op,'nDraws')
    if ismember(op.UseDist,{'PriorDraws','PostDraws'})
        op.nDraws = 1000;
    else
        op.nDraws = 1;
    end
end
if ~isfield(op,'Bands2Show'), op.Bands2Show = [50,70,90]; end
if ~isfield(op,'nPanelMax'), op.nPanelMax = 6; end
if ~isfield(op,'Vars2Show')
    list = {'ObsVar','StateVar','AuxVar'};
    jP = 0;
    for jlist=1:length(list)
        listj = list{jlist};
        for j=(1:ceil(s.n.(listj)/op.nPanelMax))
            jP = jP+1;
            op.Vars2Show(jP).Title = sprintf('%s%02.0f',listj,j);
            op.Vars2Show(jP).Names = s.(listj)(...
                (j-1)*op.nPanelMax+1:min(j*op.nPanelMax,s.n.(listj)));
        end
    end
end
if ~isfield(op,'Shocks2Show'), op.Shocks2Show = s.ShockVar; end

if ~isfield(op,'nSteps'), op.nSteps = 25; end
if ~isfield(op,'TickStep'), op.TickStep = 4; end

if ~isfield(op,'yPerSlack'), op.yPerSlack = 0.05; end
if ~isfield(op,'yMaxSlack'), op.yMaxSlack = []; end
if ~isfield(op,'yMinScale'), op.yMinScale = 0; end

if ~isfield(op,'ShowFig'), op.ShowFig = 0; end

if ~isfield(s,'PlotDir') || ~isfield(s.PlotDir,'IRF')
    s.PlotDir.IRF = 'Plots_IRF/';
end
if ~isdir(s.PlotDir.IRF), mkdir(s.PlotDir.IRF), end
if ~isfield(s.FileName,'PlotsIRF'), 
    s.FileName.PlotsIRF = sprintf('%sIRF',s.Spec); 
end

% Check options
nPanels = length(op.Vars2Show);
for jP=1:nPanels
    if ~isfield(op.Vars2Show(jP),'NamesPretty')
        op.Vars2Show(jP).NamesPretty = op.Vars2Show(jP).Names;
    end
    if ~isfield(op.Vars2Show(jP),'Scale')
        op.Vars2Show(jP).Scale = 1;
    end
end
nShocks2Show = length(op.Shocks2Show);
if ~isfield(op,'ShockSize'), op.ShockSize = ones(1,nShocks2Show); end

% Save options
s.Options.IRF = op;

%% -------------------------------------------------------------------

%% Prepare Draws
if strcmp(op.UseDist,'PriorDraws')
    xd = feval(s.FileName.PriorDraw,op.nDraws);
elseif strcmp(op.UseDist,'PostDraws')
    load(FileName.MCMCDrawsRedux,'xd')
    xd = xd(:,1:op.nDraws);
else
    if ~isfield(s.Params,op.UseDist)
        fprintf(2,'Did not recognize distribution to use. Cannot proceed.\n');
        return
    end
    xd = s.Param.(op.UseDist);
end



%% -------------------------------------------------------------------

%% Finish up
s.Status.(Action) = 1;
s.TimeElapsed.(Action) = toc-s.TimeElapsed.(Action);

%% -------------------------------------------------------------------



%% Auxilliary Variables
nPanels = length(Panels);
if nPanels>1
    TempVars = union(Vars2Show{:});
else
    TempVars = Vars2Show{1};
end
if any(ismember(TempVars,ObsVar))
    isObsMats=2;
else
    isObsMats=0;
end

%% ------------------------------------------------------------------------


%% Prepare Draws
if strcmp(UseDist,'priordraws')
    if ~isfield(FileName,'GenPriorDraw')
        MakeGenPriorDraw
    end
    xd = feval(FileName.GenPriorDraw,nDrawsUsed);
elseif strcmp(UseDist,'postdraws')
    if ~isfield(FileName,'MCMCDrawsRedux')
        MakeMCMCDrawsRedux
    end
    load(FileName.MCMCDrawsRedux,'xd')
    xd = xd(:,1:nDrawsUsed);
else
    if ~isfield(Params,UseDist)
        fprintf(2,'Did not recognize distribution to use. Cannot proceed.\n');
        return
    end
    xd = [Params(:).(UseDist)]';
end
nDrawsUsed = size(xd,2);

%% Create Indices
nShocks2Show = length(Shocks2Show);
idxShocks = zeros(nShocks2Show,1);
for ii = 1:nShocks2Show
    idxShocks(ii) = find(ismember(ShockVar,Shocks2Show(ii)));
end    
nTempVars = length(TempVars);
idxVars = zeros(nTempVars,1);
AllVars = {StateVar{:},ObsVar{:}};
for ii = 1:nTempVars
    idxVars(ii) = find(ismember(AllVars,TempVars(ii)));
end

%% Generate IRF
FileNameMats = FileName.Mats;
IRF = zeros(nTempVars,nSteps,nShocks2Show,nDrawsUsed);
fprintf('Generating IRFs...\n');
parfor j=1:nDrawsUsed
  IRF(:,:,:,j) = MakeIRFFcn(xd(:,j),nSteps,idxShocks,idxVars,FileNameMats,...
    isObsMats,nObsVar);
end

%% Reassign IRF to Vars2Show
PlotIRF = cell(1,nPanels);
for ii = 1:nPanels
    [tf,idx] = ismember(Vars2Show{ii},TempVars);
    PlotIRF{ii} = IRF(idx,:,:,:);
end

%% IRF plot prep
XTicks = 0:TickStepIRF:(nSteps-1);
for ii = 1:length(XTicks)
    XTickLabels{ii} = sprintf('%i',XTicks(ii));
end
if exist('Shape2Plot','var')
    if length(Shape2Plot)~=length(Vars2Show)
        error(...
            'Shape2Plot dimensions must be specified for each Vars2Show cell.');
    end
elseif ~exist('Shape2Plot','var')
    for ii = 1:nPanels
        Shape2Plot{ii} = [1,1];
        jDim = 1;
        while prod(Shape2Plot{ii})<length(Vars2Show{ii})
            Shape2Plot{ii}(jDim) = Shape2Plot{ii}(jDim)+1;
            jDim = ~(jDim-1)+1;
        end
    end
end

%% Plot IRFs
fprintf('Plotting IRFs...\n');
tid = 0:(nSteps-1);
for js=1:length(Shocks2Show)
    for ii = 1:nPanels
        if ShowFig 
            figure('Name',['Responses to innovation in ',Shocks2Show{js}])
        else
            figure('Visible','off')
        end
        for jj=1:length(Vars2Show{ii})
            subplot(Shape2Plot{ii}(1),Shape2Plot{ii}(2),jj);
            if nDrawsUsed==1
                vcPlot(tid,Vars2ShowScale{ii}(jj)*PlotIRF{ii}(jj,:,js))
            else
                vcPlotDistBands(tid,squeeze(...
                    Vars2ShowScale{ii}(jj)*PlotIRF{ii}(jj,:,js,:))',...
                                'Bands2Show',Bands2Show);
            end
            hold on
            plot(tid,zeros(1,nSteps),':k')
            hold off
            title(Vars2ShowPretty{ii}{jj})
            set(gca,'Xtick',XTicks,'XTickLabel',XTickLabels);
            axis tight
            yl = get(gca,'YLim');
            ySlack = max([yPerSlack*(yl(2)-yl(1)),yMaxSlack]);
            set(gca,'YLim',[min(yl(1)-ySlack,-yMinScale),...
                            max(yMinScale,yl(2)+ySlack)])
        end
        vcPrintPDF(sprintf('%s%s_%s_%s',...
                           PlotDir.IRF,FileName.PlotsIRF,...
                           Panels{ii},Shocks2Show{js}))
    end
end

%% close figures
if ~ShowFig, close all, end

