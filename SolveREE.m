function [REE,fmat,fwt,ywt,gev]=SolveREE(StateEq,Author,varargin)

% SolveREE
%
% Uses gensys to solve for the REE
%
% Option: Author (string)
% If set to 'CS' (default) it uses the original verion. If set to 'JW' it uses
% the fast gensys from Jae Won and if it does not yield a normal solution it
% runs the original gensys.
%
% ...........................................................................
% 
% Created: January 25, 2016 by Vasco Curdia
% 
% Copyright 2016 by Vasco Curdia

%% -------------------------------------------------------------------

%% Preamble
REE.GBar = [];
REE.G1 = [];
REE.G2 = [];
REE.eu = [0;0];
if ~exist('Author','var') || isempty(Author), Author = 'CS'; end

%% Run JW
if strcmp(Author,'JW')
    [REE.G1,REE.GBar,REE.G2,fmat,fwt,ywt,gev,REE.eu] = ...
        fastgensysJaeWonvb(StateEq.Gamma0,StateEq.Gamma1,StateEq.GammaBar,...
                           StateEq.Gamma2,StateEq.Gamma3,varargin{:});
end

%% Run CS
if all(REE.eu(:)==1), return, end
[REE.G1,REE.GBar,REE.G2,fmat,fwt,ywt,gev,REE.eu] = ...
    gensysvb(StateEq.Gamma0,StateEq.Gamma1,StateEq.GammaBar,...
             StateEq.Gamma2,StateEq.Gamma3,varargin{:});

%% -------------------------------------------------------------------

