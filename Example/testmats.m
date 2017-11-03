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
NSEq = {'c*d-alpha'};
CParam = {'d','4+a+b'};
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
f.NSEq = @(x)buildmat(ftmp,x,nNSParam,1);
clear ftmp

function y = solveNSParam(x)
    

function f=NumSolveEq(x)
    for jx=1:size(x,2)
        rA = x(1,jx);
        rB = x(2,jx);
        EvalCompoundParam
        f(1,jx) = (rA+rB)/2-r;
        f(2,jx) = rA+0.5/400-rB;
    end
end
NumSolveGuess = [...
    1.0000000000000000;
    1.0000000000000000;
    ];
NumSolveOptions = optimoptions(@fsolve);
NumSolveOptions.Display = 'off';
[NumSolveSolution,NumSolveResidual,NumSolveRC,NumSolveOutput] = fsolve(@NumSolveEq,NumSolveGuess,NumSolveOptions);


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
