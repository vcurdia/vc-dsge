function maxpostreport(obj,varargin)

% maxpostreport
% 
% make report on maxpost results
%
% see also:
% DSGE.Model, maxpost
%
% ...
%
% Created: November 20, 2024
% Copyright (c) 2024 Vasco Curdia

    %% preamble
    ReportFileName = sprintf('%s-param-postmode',obj.Name);
    ReportTitle = sprintf('%s\\\\[30pt]Parameter Analysis\\\\Posterior Mode',obj.Name);

    %% Options
    op.Table = DSGE.Options.Table;
    op = updateoptions(op,varargin{:});

    mats = obj.mats(obj.Post.Mode,'SolveREE',0);
    xAux = mats.AuxParam;


    %% show results on screen
    pNames = obj.Param.Names;
    pNameLength = [cellfun('length',pNames)];
    nameLengthMax = max(pNameLength);
    auxN = obj.AuxParam.N;
    auxNames = obj.AuxParam.Names;
    auxNameLength = [cellfun('length',auxNames)];
    nameLengthMax = max(nameLengthMax,max(auxNameLength));

    fprintf('\nResults from maximization of posterior:')
    fprintf('\n=======================================\n')
    DispList = {'','',pNames;
                'Prior','Dist',obj.Prior.Dist;
                '','   Mode',obj.Prior.Mode;
                '','   Mean',obj.Prior.Mean;
                '','     SD',obj.Prior.SD;
                '','     5%',obj.Prior.Prc05;
                '',' Median',obj.Prior.Median;
                '','    95%',obj.Prior.Prc95;
                'Posterior','   Mode',obj.Post.Mode;
                '','     SD',obj.Post.SD;
               };
    nc = size(DispList,1);
    for jr=1:2
        str2show = sprintf(['%-',int2str(nameLengthMax),'s'],DispList{1,jr});
        str2show = sprintf('%s  %-5s',str2show,DispList{2,jr});
        for jc=3:nc
            str2show = sprintf('%s  %-7s',str2show,DispList{jc,jr});
        end
        disp(str2show)
    end
    for jp=1:obj.Param.N
        str2show = sprintf(['%',int2str(nameLengthMax),'s'],DispList{1,3}{jp});
        str2show = sprintf('%s  %5s',str2show,DispList{2,3}{jp});
        for jc=3:nc
            str2show = sprintf('%s  %7.3f',str2show,DispList{jc,3}(jp));
        end
        disp(str2show)
    end
    fprintf('\n')

    DispList = {'','',auxNames;
                'Prior','   Mean',obj.Prior.Sample.AuxParam.Mean;
                '','     5%',obj.Prior.Sample.AuxParam.Prc05;
                '',' Median',obj.Prior.Sample.AuxParam.Median;
                '','    95%',obj.Prior.Sample.AuxParam.Prc95;
                'Posterior','   Mode',xAux;
               };
    nc = size(DispList,1);
    for jr=1:2
        str2show = sprintf(['%-',int2str(nameLengthMax),'s'],DispList{1,jr});
        for jc=2:nc
            str2show = sprintf('%s  %-7s',str2show,DispList{jc,jr});
        end
        disp(str2show)
    end
    for jp=1:auxN
        str2show = sprintf(['%',int2str(nameLengthMax),'s'],DispList{1,3}{jp});
        for jc=2:nc
            str2show = sprintf('%s  %7.3f',str2show,DispList{jc,3}(jp));
        end
        disp(str2show)
    end
    fprintf('\n')

    fprintf('\nposterior log-pdf at mode: %.6f\n\n',obj.Post.ModeLPDF)

    %% create report
    fprintf('Making report: %s\n',ReportFileName);
    fid = createtex(ReportFileName,ReportTitle);

    fprintf(fid,'\\begin{equation*} \n');
    fprintf(fid,'\\begin{tabular}{rl} \n');
    fprintf(fid,'posterior log density at mode: & %.4f\n',obj.Post.ModeLPDF);
    fprintf(fid,'\\end{tabular}\n');
    fprintf(fid,'\\end{equation*}\n');

    fprintf(fid,'\\newpage \n');
    np = obj.Param.N;
    str = [' & $%.',int2str(op.Table.Precision),'f$'];
    tableBreaks = settablebreaks(np,op.Table.MaxRows);
    idxPar = 0;
    nBreaks = length(tableBreaks);
    for jBreak=1:nBreaks
        idxPar = (idxPar(end)+1):tableBreaks(jBreak);
        if nBreaks==1
            fprintf(fid,'\\subsection{Parameters}\n');
        else
            fprintf(fid,'\\subsection{Parameters (%.0f/%.0f)}\n',...
                    jBreak,nBreaks);
        end
        fprintf(fid,'\\begin{equation*}\n');
        if op.Table.MoveLeft
            fprintf(fid,'\\hspace{-0.5in}\n');
        end
        fprintf(fid,'\\begin{tabular}{l%s} \n',repmat('r',1,1+7+1+2));
        fprintf(fid,'\\hline\\hline\\\\[-1.5ex]\n');
        fprintf(fid,'& \\multicolumn{7}{c}{Prior} ');
        fprintf(fid,'& & \\multicolumn{2}{c}{Posterior} \\\\[0.5ex]\n');
        fprintf(fid,'& Dist & Mode & Mean & SD & 5\\%% & Median & 95\\%% ');
        fprintf(fid,'& & Mode & SD \n');
        fprintf(fid,'\\\\[0.5ex]\\hline\\\\[-1.5ex]\n');
        for jr=idxPar
            fprintf(fid,'%s',obj.Param.PrettyNames{jr});
            fprintf(fid,' & %s', obj.Prior.Dist{jr});
            fprintf(fid,str,obj.Prior.Mode(jr));
            fprintf(fid,str,obj.Prior.Mean(jr));
            fprintf(fid,str,obj.Prior.SD(jr));
            fprintf(fid,str,obj.Prior.Prc05(jr));
            fprintf(fid,str,obj.Prior.Median(jr));
            fprintf(fid,str,obj.Prior.Prc95(jr));
            fprintf(fid,' &');
            fprintf(fid,str,obj.Post.Mode(jr));
            fprintf(fid,str,obj.Post.SD(jr));
            fprintf(fid,' \\\\\n');
            if ismember(jr,op.Table.Lines) && jr~=idxPar(end)
                fprintf(fid,'\\\\[-1.5ex]\\hline\\\\[-1.5ex]\n');
            end        
        end
        fprintf(fid,'\\\\[-1.5ex]\\hline\\hline\n');
        fprintf(fid,'\\end{tabular}\n');
        fprintf(fid,'\\end{equation*}\n');
        fprintf(fid,'\\clearpage\n');
    end

    tableBreaks = settablebreaks(auxN,op.Table.MaxRows);
    idxPar = 0;
    nBreaks = length(tableBreaks);
    for jBreak=1:nBreaks
        idxPar = (idxPar(end)+1):tableBreaks(jBreak);
        if nBreaks==1
            fprintf(fid,'\\subsection{Auxiliary Parameters}\n');
        else
            fprintf(fid,'\\subsection{Auxiliary Parameters (%.0f/%.0f)}\n',...
                    jBreak,nBreaks);
        end
        fprintf(fid,'\\begin{equation*}\n');
        fprintf(fid,'\\begin{tabular}{l%s} \n',repmat('r',1,1+4+1+1));
        fprintf(fid,'\\hline\\hline\\\\[-1.5ex]\n');
        fprintf(fid,'& \\multicolumn{4}{c}{Prior} & & Posterior\\\\[0.5ex]\n');
        fprintf(fid,'& Mean & 5\\%% & Median & 95\\%% & & Mode \n');
        fprintf(fid,'\\\\[0.5ex]\\hline\\\\[-1.5ex]\n');
        for jr=idxPar
            fprintf(fid,'%s',obj.AuxParam.PrettyNames{jr});
            fprintf(fid,str,obj.Prior.Sample.AuxParam.Mean(jr));
            fprintf(fid,str,obj.Prior.Sample.AuxParam.Prc05(jr));
            fprintf(fid,str,obj.Prior.Sample.AuxParam.Median(jr));
            fprintf(fid,str,obj.Prior.Sample.AuxParam.Prc95(jr));
            fprintf(fid,' &');
            fprintf(fid,str,xAux(jr));
            fprintf(fid,' \\\\\n');
            if ismember(jr,op.Table.Lines) && jr~=idxPar(end)
                fprintf(fid,'\\\\[-1.5ex]\\hline\\\\[-1.5ex]\n');
            end        
        end
        fprintf(fid,'\\\\[-1.5ex]\\hline\\hline\n');
        fprintf(fid,'\\end{tabular}\n');
        fprintf(fid,'\\end{equation*}\n');
        fprintf(fid,'\\clearpage\n');
    end

    fprintf(fid,'\\end{document}\n');
    fclose(fid);
    pdflatex(ReportFileName)

end

