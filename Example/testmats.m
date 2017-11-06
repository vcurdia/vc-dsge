% Test mats in object

%% Preamble
clear all
setpath
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
tt = TimeTracker;
tt.start('Setup')


%% simple example setup
Param = {'a';'b';'c'};
NSParam = {'alpha'};
NSEq = {'a-alpha'};
CParam = {'d','4+alpha+b'};
Var = {'x';'y';'z'};
Eq = {...
    'a*x-z+alpha';
    '2+x*d+b*y';
    '5-c*y';
    };
x0 = [1;2;3];

%% generate mats
AuxParam = [NSParam;CParam(:,1)];
AllParam = [Param;AuxParam];
vcsym(AllParam{:})

nNSParam = length(NSParam);
nCParam = size(CParam,1);
nAuxParam = length(AuxParam);

ftmp = cell(nCParam,1);
for j=1:nCParam
    fj = eval(CParam{j,2});
    ftmp{j} = matlabFunction(fj,'Vars',[Param;NSParam]);
end
f.CParam = @(x)buildmat(ftmp,x,nCParam,1);
clear ftmp

ftmp = cell(nNSParam,1);
for j=1:nNSParam
    fj = eval(NSEq{j});
    ftmp{j} = matlabFunction(fj,'Vars',AllParam);
end
f.NSEq = @(x)buildmat(ftmp,[x;f.CParam(x)],nNSParam,1);
clear ftmp

f.AuxParam = @(x)fAuxParam(x,f,nNSParam);

f.AllParam = @(x)[x;f.AuxParam(x)];

vcsym(Var{:})
nVar = length(Var);
var = sym(zeros(nVar,1));
for j=1:nVar
    var(j) = eval(Var{j});
end
nEq = length(Eq);
eq = sym(zeros(nEq,1));
for j=1:nEq
    eq(j) = eval(Eq{j});
end
EqMat = jacobian(eq,var);
fEq = cell(nEq*nVar,1);
for j=1:nEq*nVar
    fEq{j} = matlabFunction(EqMat(j),'Vars',AllParam);
end
f.Eq = @(x)buildmat(fEq,x,nEq,nVar);
clear fEq

%% Test functions
mats.Param = x0;
mats.AuxParam = f.AuxParam(x0);
mats.AllParam = [mats.Param;mats.AuxParam];
mats.Eq = f.Eq(mats.AllParam);

%% Finish up
tt.stop('Setup')
save testmatsws

function xaux=fAuxParam(x,f,nNSParam,xns0)
    for jx=1:size(x,2)
%         eqj = @(xns)f.NSEq([x(:,j);xns]);
        if nargin<4 || isempty(xns0)
            xns0 = ones(nNSParam,1);
        end
        NumSolveOptions = optimoptions(@fsolve);
        NumSolveOptions.Display = 'off';
        [xns1,NumSolveResidual,NumSolveRC,NumSolveOutput] = ...
            fsolve(@(xns)f.NSEq([x(:,jx);xns]),xns0,NumSolveOptions);
        xaux(:,jx) = [xns1;f.CParam([x(:,jx);xns1])];
    end
end
