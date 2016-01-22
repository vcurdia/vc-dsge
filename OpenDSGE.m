function s = OpenDSGE(Spec)

%% -------------------------------------------------------------------

if ~isdir(Spec)
    mkdir(Spec)
end
cd(Spec)

if exist([Spec,'.mat'],'file')
    fprintf('\nOpening DSGE: %s\n',Spec);
    s = load(Spec);
else
    fprintf('\nInitiating DSGE: %s\n',Spec);
    s.Spec = Spec;
    save(Spec,'-struct','s')
end

