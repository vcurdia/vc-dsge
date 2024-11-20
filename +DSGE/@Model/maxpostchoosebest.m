function maxpostchoosebest(obj,idxMax)

% maxpostchoosebest
% 
% select the maxpost run with best mode candidate
%
% see also:
% DSGE.Model, maxpost
%
% ...
%
% Created: November 20, 2024
% Copyright (c) 2024 Vasco Curdia

    pIdx = obj.Post.EstimateIdx;

    %% load minimizations output
    load([obj.Name,'-maxpost-out'],'MaxPostOut')

    %% extract the best one
    if nargin<2 || isempty(idxMax)
        [LPDFMode, idxMax] = min([MaxPostOut(:).f]);
    end
    fprintf('Using minimization %.0f\n',idxMax)

    %% update Post
    obj.Post.Mode(pIdx) = MaxPostOut(idxMax).x;
    obj.Post.ModeLPDF = -MaxPostOut(idxMax).f;
    obj.Post.Var(pIdx,pIdx) = MaxPostOut(idxMax).H;
    obj.Post.SD = diag(obj.Post.Var).^(1/2);

    %% save back minimizations with updated best choice
    save([obj.Name,'-maxpost-out'],'MaxPostOut','idxMax')

end

