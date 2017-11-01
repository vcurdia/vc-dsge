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
AuxParam = {'d','a+b'};

x0 = [1;2;3];

%% generate mats
AllParam = [Param;AuxParam(:,1)];
vcsym(AllParam{:})

mats.Param = @(x)x;

nAuxParam = size(AuxParam,1);
f.AuxParam = cell(nAuxParam,1);
for j=1:nAuxParam
    Auxj = eval(AuxParam{j,2});
    f.AuxParam{j} = matlabFunction(Auxj,'Vars',Param);
end
mats.AuxParam = @(x)buildmat(f.AuxParam,x,nAuxParam,1);


%% Finish up
tt.stop('Setup')
