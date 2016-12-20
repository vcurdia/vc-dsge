function obj = GenPost(obj)

% GenPost
%
% Generates the posterior function for the model
%
% See also:
% DSGE, SetupMyDSGE
%
% ...........................................................................
%
% Created: December 19, 2016 by Vasco Curdia
% 
% Copyright (C) 2016 Vasco Curdia

%% -------------------------------------------------------------------

%% Preamble

action = 'GenPost';
obj = obj.TrackTime(action,1);

fprintf('\n*** Generate posterior function\n')

%% -------------------------------------------------------------------

%% Initiate file

obj.FileName.Post = sprintf('%s_Post',obj.Name);

fid = fopen([obj.FileName.Post,'.m'],'wt');
fprintf(fid,'function post = %s(x,data,varargin)\n\n',obj.FileName.Post);
fprintf(fid,'%% Created: %.0f/%.0f/%.0f %.0f:%.0f:%.0fs\n',clock);

%% -------------------------------------------------------------------

fprintf(fid,'\n%% Default options\n');
fprintf(fid,'op.fid = 1;\n');
fprintf(fid,'op.verbose = 0;\n');

fprintf(fid,'\n%% Update options\n');
fprintf(fid,'if length(varargin)>0 && isstruct(varargin{1})\n');
fprintf(fid,'    opnew = varargin{1};\n');
fprintf(fid,'    opnewfields = fieldnames(opnew);\n');
fprintf(fid,'    for jop=1:length(opnewfields)\n');
fprintf(fid,'        op.(opnewfields{jop}) = opnew.(opnewfields{jop});\n');
fprintf(fid,'    end\n');
fprintf(fid,'    varargin(1) = [];\n');
fprintf(fid,'end\n');
fprintf(fid,'for jop=1:(length(varargin)/2)\n'); 
fprintf(fid,'    op.(varargin{(jop-1)*2+1}) = varargin{jop*2};\n');
fprintf(fid,'end\n');

fprintf(fid,'\n%% Map parameters\n\n');
for j=1:obj.Param.N
    fprintf(fid,'%s = x(%.0f);\n',obj.Param.Names{j},j);
end

fprintf(fid,'\n%% Map fixed parameters\n\n');
fprintf(fid,'xFix = [...\n');
fprintf(fid,'    %.16f;\n',obj.FixParam.Values);
fprintf(fid,'    ];\n');

fprintf(fid,'\n%% Evaluate priors\n\n');
fprintf(fid,'post = 0;\n');
for j=1:obj.Param.N
    fprintf(fid,'post = post + %s;\n',obj.Param.PriorLPdfCmd{j});
end
fprintf(fid,'post = post + %.16f;\n',obj.Prior.LogTruncationCorrection);
fprintf(fid,'if post==-inf, post = inf; return, end;\n\n');

fprintf(fid,'\n%% Get Mats\n');
fprintf(fid,'Mats = %s([x;xFix],op);\n',obj.FileName.Mats);
fprintf(fid,'if ~all(Mats.REE.eu==1) || Mats.KF.sig00rc~=0\n');
fprintf(fid,'    post = inf;\n');
fprintf(fid,'    return\n');
fprintf(fid,'end\n\n');

fprintf(fid,'\n%% Kalman Filter\n\n');
fprintf(fid,'stt = Mats.KF.s00;\n');
fprintf(fid,'sigtt = Mats.KF.sig00;\n');
fprintf(fid,'StateVartt = zeros(%.0f,%.0f);\n',obj.StateVar.N,obj.T);
fprintf(fid,'SIGtt = zeros(%.0f,%.0f,%.0f);\n',obj.StateVar.N,obj.StateVar.N,...
        obj.T);
fprintf(fid,'DataDetrended = data-Mats.KF.ObsVarBar'';\n');
fprintf(fid,'for t=1:%.0f;\n',obj.T);
fprintf(fid,'    idxNoNaN = ~isnan(data(t,:));\n');
fprintf(fid,'    [stt,sigtt,lh,ObsVarhat]=kf(...\n');
fprintf(fid,'        DataDetrended(t,idxNoNaN)'',...\n');
fprintf(fid,'        Mats.ObsEq.H(idxNoNaN,:),stt,sigtt,...\n');
fprintf(fid,'        Mats.REE.G1,Mats.REE.G2);\n');
fprintf(fid,'    if t>%.0f\n',obj.NPreSample);
fprintf(fid,'        post = post + lh*[1;1];\n');
fprintf(fid,'    end\n');
fprintf(fid,'    StateVartt(:,t) = stt;\n');
fprintf(fid,'    SIGtt(:,:,t) = sigtt;\n');
fprintf(fid,'end\n\n');

fprintf(fid,'%% Add normalization\n');
fprintf(fid,['post = -( post ',...
             '- sum(sum(~isnan(data(%.0f:%.0f,:))))/2*log(2*pi) );\n\n'],...
             obj.NPreSample+1,obj.T);

%% -------------------------------------------------------------------

%% close file
fprintf(fid,'end\n');
fclose(fid);
    
%% -------------------------------------------------------------------

%% Test Posterior
fprintf('Testing posterior function...\n');
post = -feval(obj.FileName.Post,obj.Param.PriorMean,obj.Data);
%post = -feval(obj.FileName.Post,obj.Prior.UnconstrainedParam.PriorMean,obj.Data);
fprintf('The log-posterior value using the prior mean is %0.4f.\n',post);

%% -------------------------------------------------------------------

%% Finish up
obj = obj.TrackTime(action,0);


%% -------------------------------------------------------------------
