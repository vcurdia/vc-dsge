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
AuxParam = {'d','4+a+b'};
Var = {'x';'y';'z'};
Eq = {...
    'a*x-z';
    '2+x*d+b*y';
    '5-c*y';
    };
x0 = [1;2;3];

%% generate mats
AllParam = [Param;AuxParam(:,1)];
vcsym(AllParam{:})

nAuxParam = size(AuxParam,1);
fAuxParam = cell(nAuxParam,1);
for j=1:nAuxParam
    Auxj = eval(AuxParam{j,2});
    fAuxParam{j} = matlabFunction(Auxj,'Vars',Param);
end
f.AuxParam = @(x)buildmat(fAuxParam,x,nAuxParam,1);
clear fAuxParam

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
