function OpenDSGE(Spec)

%% -------------------------------------------------------------------

if ~isdir(Spec)
    mkdir(Spec)
end
cd(Spec)

if exist('Workspace','file')
    fprintf('\nOpening DSGE: %s\n',Spec);
    load('Workspace')
else
    fprintf('\nInitiating DSGE: %s\n',Spec);
end

