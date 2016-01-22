function s = CloseDSGE(s,BaseFolder)

fprintf('\nSaving workspace and exiting DSGE folder\n')
save(s.Spec,'-struct','s')
if nargin>1
    cd(BaseFolder)
else
    cd ..
end

