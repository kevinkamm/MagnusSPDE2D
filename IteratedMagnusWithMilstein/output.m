function output(fileName,methodName,...
                xl,xu,nx,vl,vu,nv,T,tnEuler,tnMagnus,dtEulerRef,dtEuler,dtMagnus,M,...
                h,fx,fv,gxx,gxv,gvv,sx,sv,phi,...
                ctimeExactTotal,ctimeEulerRefTotal,ctimeEulerTotal,ctimeMagnus1Total,ctimeMagnus2Total,ctimeMagnus3Total,...
                xgrid,vgrid,kappas,...
                errDisExactEulerRef,errDisExactEuler,errDisExactMagnus1,errDisExactMagnus2,errDisExactMagnus3,...
                errDisEulerRefEuler,errDisEulerRefMagnus1,errDisEulerRefMagnus2,errDisEulerRefMagnus3,...
                trajectoryPlots,cdfExactPlots,cdfEulerRefPlots,patterns,...
                hasExact,useExact,useEulerRef,useEuler,useMagnus1,useMagnus2,useMagnus3,...
                tEval,absErrorExactVid,absErrorEulerRefVid,errorDisExactVid,errorDisEulerRefVid,...
                errAbsAvgExactEulerRef,errAbsAvgExactEuler,errAbsAvgExactMagnus1,errAbsAvgExactMagnus2,errAbsAvgExactMagnus3,...
                errAbsAvgEulerRefEuler,errAbsAvgEulerRefMagnus1,errAbsAvgEulerRefMagnus2,errAbsAvgEulerRefMagnus3,...
                varargin)
backgroundColor='w';
compileLatex = true;
for iV=1:2:length(varargin)
    switch varargin{iV}
        case 'backgroundColor'
            backgroundColor = varargin{iV+1};
        case 'compileLatex'
            compileLatex = varargin{iV+1};
    end
end
fclose('all');
picType='eps';
saveParam='epsc';
numDict={'One','Two','Three','Four','Five','Six','Seven','Eight','Nine'};

root=[pwd, '/' ,'Results'];
pdfRoot=[root,'/','Pdf'];
tempPath=[pdfRoot,'/','temp'];
copyPath=[pdfRoot,'/',methodName,'/',fileName];
templatePath=[tempPath,'/','template', '.','tex'];
if compileLatex
    outputFilePath={[copyPath,'/','template','.','pdf'],...
                [copyPath,'/','template','.','tex']};
    copyFilePath={[copyPath,'/',fileName,'.','pdf'],...
                [copyPath,'/',fileName,'.','tex']};
else
    outputFilePath={[copyPath,'/','template','.','tex']};
    copyFilePath={[copyPath,'/',fileName,'.','tex']};
end
inputPath=[tempPath,'/','input','.','tex'];

% Delete auxiliary files in temp folder
try
cleanDir(tempPath,{'template.tex'});
catch
    disp('Error in cleandir')
end
delDir(copyPath)
mkDir(copyPath)

inputFile=fopen(inputPath,'w+');
%% Head
fprintf(inputFile,...
        '\\section{2D Magnus}\n');
fprintf(inputFile,...
    ['\\begin{align*}\n\\left\\{\n\\begin{aligned}[c]\\arraycolsep=0pt\n&\n\\begin{aligned}[t]\\arraycolsep=0pt\n',...
    '&\\biggl(\n-\\partial_t u_t(x,v) \n+ h(x,v) u_t(x,v)\n+ f^x(x,v) \\partial_x u_t(x,v)\n',...
    '+ f^v(x,v) \\partial_v u_t(x,v)\n\\\\&\\quad\n+ \\frac{1}{2} g^{xx}(x,v) \\partial_{xx} u_t(x,v)\n',...
    '+ g^{xv}(x,v) \\partial_{xv} u_t(x,v)\n+ \\frac{1}{2} g^{vv}(x,v) \\partial_{vv} u_t(x,v)\n',...
    '\\biggr) dt\n\\\\&\\quad\n+ \\left(\\sigma^{x}(x,v) \\partial_{x} u_t(x,v) \n',...
    '+ \\sigma^{v}(x,v) \\partial_{v} u_t(x,v)\\right) dW_t\n= 0\n\\end{aligned}\\\\&\n',...
    'u_0(x,v)=\\phi(x,v)\n\\end{aligned}\n\\right.\n\\end{align*}\n']);
%% Model Coefficient Functions
fprintf(inputFile,...
        '\\subsection{Model Coefficient Functions}\n');
[latexFilePath,latexCommand]=saveModelCoefficients('ModelCoeff');
fprintf(inputFile,...
        '\t\\input{%s}\n',changeSlash(latexFilePath));
for iLC=1:1:length(latexCommand)
    fprintf(inputFile,...
            '\t%s\n\n\n',latexCommand{iLC});
end
%% Model Paramters
fprintf(inputFile,...
        '\\subsection{Numerical Parameters}\n');
[latexFilePath,latexCommand]=saveNumericalParameters('NumericParam');
fprintf(inputFile,...
        '\t\\input{%s}\n',changeSlash(latexFilePath));
for iLC=1:1:length(latexCommand)
    fprintf(inputFile,...
            '\t%s\n\n\n',latexCommand{iLC});
end
%% Computational Times
fprintf(inputFile,...
        '\\subsection{Computational Times}\n');
[latexFilePath,latexCommand]=saveComputationalTimes('CompTimes');
fprintf(inputFile,...
        '\t\\input{%s}\n',changeSlash(latexFilePath));
for iLC=1:1:length(latexCommand)
    fprintf(inputFile,...
            '\t%s\n\n\n',latexCommand{iLC});
end
%% Errors
fprintf(inputFile,...
        '\\subsection{Errors}\n');
% Error regions
fprintf(inputFile,...
        '\\subsubsection{Error regions}\n');
[latexFilePath,latexCommand]=saveErrorRegion('Errors');
fprintf(inputFile,...
        '\t\\input{%s}\n',changeSlash(latexFilePath));
for iLC=1:1:length(latexCommand)
    fprintf(inputFile,...
            '\t%s\n\n\n',latexCommand{iLC});
end
% Absolute errors
fprintf(inputFile,...
        '\\subsubsection{Abs Errors}\n');
[latexFilePath,latexCommand]=saveAbsErrors('AbsErrors');
fprintf(inputFile,...
        '\t\\input{%s}\n',changeSlash(latexFilePath));
for iLC=1:1:length(latexCommand)
    fprintf(inputFile,...
            '\t%s\n\n\n',latexCommand{iLC});
end
% Error distributions
fprintf(inputFile,...
        '\\subsubsection{Error Distributions}\n');
[latexFilePath,latexCommand]=saveErrorDis('ErrorDis');
fprintf(inputFile,...
        '\t\\input{%s}\n',changeSlash(latexFilePath));
for iLC=1:1:length(latexCommand)
    fprintf(inputFile,...
            '\t%s\n\n\n',latexCommand{iLC});
end
%% Error videos
fprintf(inputFile,...
        '\\subsection{Plots of Mean Absolute Errors}\n');
if hasExact && useExact && ~isempty(absErrorExactVid)
    fprintf(inputFile,...
        '\\subsection{Reference method: exact}\n');
    latexFilePath=saveVideos(absErrorExactVid,'absErrorExactVideo');
    fprintf(inputFile,...
            '\t\\input{%s}\n',changeSlash(latexFilePath));
end
if useEulerRef && ~isempty(absErrorEulerRefVid)
    fprintf(inputFile,...
        '\\subsection{Reference method: Euler Ref}\n');
    latexFilePath=saveVideos(absErrorEulerRefVid,'absErrorEulerRefVideo');
    fprintf(inputFile,...
            '\t\\input{%s}\n',changeSlash(latexFilePath));
end
if hasExact && useExact && ~isempty(errorDisExactVid)
%     fprintf(inputFile,...
%         '\\subsection{Reference method: exact}\n');
    latexFilePath=saveVideos(errorDisExactVid,'errorDisExactVid');
%     fprintf(inputFile,...
%             '\t\\input{%s}\n',changeSlash(latexFilePath));
end
if useEulerRef && ~isempty(errorDisEulerRefVid)
%     fprintf(inputFile,...
%         '\\subsection{Reference method: Euler Ref}\n');
    latexFilePath=saveVideos(errorDisEulerRefVid,'errorDisEulerRefVid');
%     fprintf(inputFile,...
%             '\t\\input{%s}\n',changeSlash(latexFilePath));
end
%% Plots of Trajectories
if ~isempty(trajectoryPlots)
    fprintf(inputFile,...
            '\\subsection{Plots of Trajectories}\n');
    latexFilePath=saveFiguresLandscape(trajectoryPlots,'TR',['TR']);
    fprintf(inputFile,...
            '\t\\input{%s}\n',changeSlash(latexFilePath));
end
%% Plots of CDFs
fprintf(inputFile,...
            '\\subsection{Plots of CDFs}\n');
if ~isempty(cdfExactPlots)
    fprintf(inputFile,...
            '\\subsubsection{Exact as reference}\n');
    latexFilePath=saveFiguresLandscape(cdfExactPlots,'CDF',['CDFExact']);
    fprintf(inputFile,...
            '\t\\input{%s}\n',changeSlash(latexFilePath));
end
if ~isempty(cdfEulerRefPlots)
    fprintf(inputFile,...
            '\\subsubsection{Euler Ref as reference}\n');
    latexFilePath=saveFiguresLandscape(cdfEulerRefPlots,'CDF',['CDFEulerRef']);
    fprintf(inputFile,...
            '\t\\input{%s}\n',changeSlash(latexFilePath));
end
%% Plots of Sparsity Patterns
    if ~isempty(patterns)
    fprintf(inputFile,...
            '\\subsection{Sparsity Patterns}\n');
    latexFilePath=saveFigures(patterns,'SparsityPatterns',['SPP']);
    fprintf(inputFile,...
            '\t\\input{%s}\n',changeSlash(latexFilePath));
    end
%% Close file
fclose(inputFile);
%% Compile Latex
if compileLatex
currFolder=cd(tempPath);
str1=sprintf('pdflatex %s',templatePath);
% system(str1)
[returncode, ~] = system(str1);
cd(currFolder);
end

copyfile(tempPath,copyPath);
% Renaming files
for iFile=1:1:length(copyFilePath)
    movefile(outputFilePath{iFile},copyFilePath{iFile})
end
% Delete auxiliary latex files in copy folder
delete([copyPath,'/','template*.*']);
% Delete auxiliary files in temp folder
% cleanDir(tempPath,{'template.tex'});
fclose('all');
%% Latex functions
    function [latexFilePath,latexCommand]=saveErrorRegion(saveAt)
            latexCommand={};
            if strcmp(saveAt, '')
                latexFilePath='regions.tex';
            else
                latexFilePath=[saveAt,'/','regions.tex'];
                mkDir([tempPath,'/',saveAt]);
            end
            file = fopen([tempPath,'/',latexFilePath],'w+');
            latexCommand{end+1}=['\errorRegion'];
            fprintf(file,'\\newcommand{%s}{\n',latexCommand{end});
            fprintf(file,'\\begin{compactenum}\n');
            for iKappa = 1:1:length(kappas)
                fprintf(file,'\\item $\\kappa=%d$: ',kappas(iKappa));
                kappa=2^iKappa; 
                regionX=1+floor(nx/2-nx/(2*kappa)):1:ceil(nx/2+nx/(2*kappa));
                regionV=1+floor(nv/2-nv/(2*kappa)):1:ceil(nv/2+nv/(2*kappa));
                fprintf(file,'$\\left[%1.3f,%1.3f\\right]\\times\\left[%1.3f,%1.3f\\right]$\n',...
                    xgrid(regionX(1)),xgrid(regionX(end)),vgrid(regionV(1)),vgrid(regionV(end)) );
            end
            fprintf(file,'\\end{compactenum}\n');
            fprintf(file,'}\n');
            fclose(file);
    end
    function [latexFilePath,latexCommand]=saveAbsErrors(saveAt)
        latexCommand={};
        if strcmp(saveAt, '')
            latexFilePath='absErrors.tex';
        else
            latexFilePath=[saveAt,'/','absErrors.tex'];
            mkDir([tempPath,'/',saveAt]);
        end
        file = fopen([tempPath,'/',latexFilePath],'w+');
        if useExact && hasExact
            latexCommand{end+1}=['\absErrorsExact'];
            fprintf(file,'\\newcommand{%s}{\n',latexCommand{end});
            for tk=1:1:length(tEval)
                if mod(tEval(tk),0.1)==0
                    fprintf(file,...
                        'Average Absolute Error at time %1.3f\n\n\n',tEval(tk));
                    corner='\\diagbox{Method}{$\\kappa$}';
                    head=kappas;
                    body=[];
                    index={};
                    if useEulerRef
                        body=cat(1,body,errAbsAvgExactEulerRef{tk}');
                        index{end+1}='Euler Ref';
                    end
                    if useEuler
                        body=cat(1,body,errAbsAvgExactEuler{tk}');
                        index{end+1}='Euler';
                    end
                    if useMagnus1
                        body=cat(1,body,errAbsAvgExactMagnus1{tk}');
                        index{end+1}='m1';
                    end
                    if useMagnus2
                        body=cat(1,body,errAbsAvgExactMagnus2{tk}');
                        index{end+1}='m2';
                    end
                    if useMagnus3
                        body=cat(1,body,errAbsAvgExactMagnus3{tk}');
                        index{end+1}='m3';
                    end
                    fprintf(file,'\\begin{tabular}{@{}c|*{%d}{c}@{}}\n',length(kappas));
                    fprintf(file,...
                        mat2Table(corner,head,index',body,{'','%1.3e'},{},...
                        'headHook','\\toprule\\\\\n'));
                    fprintf(file,'\\end{tabular}\n\n\n');
                end
            end
            fprintf(file,'}\n');
        end
        if useEulerRef && ~(useExact && hasExact)
            latexCommand{end+1}=['\absErrorsEulerRef'];
            fprintf(file,'\\newcommand{%s}{\n',latexCommand{end});
            for tk=1:1:length(tEval)
                if mod(tEval(tk),0.1)==0
                    fprintf(file,...
                        'Average Absolute Error at time %1.3f\n\n\n',tEval(tk));
                    corner='\\diagbox{Method}{$\\kappa$}';
                    head=kappas;
                    body=[];
                    index={};
                    if useEuler
                        body=cat(1,body,errAbsAvgEulerRefEuler{tk}');
                        index{end+1}='Euler';
                    end
                    if useMagnus1
                        body=cat(1,body,errAbsAvgEulerRefMagnus1{tk}');
                        index{end+1}='m1';
                    end
                    if useMagnus2
                        body=cat(1,body,errAbsAvgEulerRefMagnus2{tk}');
                        index{end+1}='m2';
                    end
                    if useMagnus3
                        body=cat(1,body,errAbsAvgEulerRefMagnus3{tk}');
                        index{end+1}='m3';
                    end
                    fprintf(file,'\\begin{tabular}{@{}c|*{%d}{c}@{}}\n',length(kappas));
                    fprintf(file,...
                        mat2Table(corner,head,index',body,{'','%1.3e'},{},...
                        'headHook','\\toprule\\\\\n'));
                    fprintf(file,'\\end{tabular}\n\n\n');
                end
            end
            fprintf(file,'}\n');
        end
        fclose(file);
    end
    function [latexFilePath,latexCommand]=saveErrorDis(saveAt)
        latexCommand={};
        if strcmp(saveAt, '')
            latexFilePath='errorDis.tex';
        else
            latexFilePath=[saveAt,'/','errorDis.tex'];
            mkDir([tempPath,'/',saveAt]);
        end
        file = fopen([tempPath,'/',latexFilePath],'w+');
        if useExact && hasExact
            latexCommand{end+1}=['\errDisExact'];
            fprintf(file,'\\newcommand{%s}{\n',latexCommand{end});
            for tk=1:1:length(tEval)
                if mod(tEval(tk),0.1)==0
                    fprintf(file,...
                        'Error distributions at time %1.3f\n\n\n',tEval(tk));
                    corner='\\diagbox{Method}{$\\kappa$}';
                    head=kappas;
                    body=[];
                    index={};
                    if useEulerRef
                        temp1=errDisExactEulerRef{tk};
                        temp=zeros(1,length(temp1));
                        for it1 = 1:1:length(temp1)
                            temp(it1)=mean(temp1{it1},1);
                        end
                        body=cat(1,body,temp);
                        index{end+1}='Euler Ref';
                    end
                    if useEuler
                        temp1=errDisExactEuler{tk};
                        temp=zeros(1,length(temp1));
                        for it1 = 1:1:length(temp1)
                            temp(it1)=mean(temp1{it1},1);
                        end
                        body=cat(1,body,temp);
                        index{end+1}='Euler';
                    end
                    if useMagnus1
                        temp1=errDisExactMagnus1{tk};
                        temp=zeros(1,length(temp1));
                        for it1 = 1:1:length(temp1)
                            temp(it1)=mean(temp1{it1},1);
                        end
                        body=cat(1,body,temp);
                        index{end+1}='m1';
                    end
                    if useMagnus2
                        temp1=errDisExactMagnus2{tk};
                        temp=zeros(1,length(temp1));
                        for it1 = 1:1:length(temp1)
                            temp(it1)=mean(temp1{it1},1);
                        end
                        body=cat(1,body,temp);
                        index{end+1}='m2';
                    end
                    if useMagnus3
                        temp1=errDisExactMagnus3{tk};
                        temp=zeros(1,length(temp1));
                        for it1 = 1:1:length(temp1)
                            temp(it1)=mean(temp1{it1},1);
                        end
                        body=cat(1,body,temp);
                        index{end+1}='m3';
                    end
                    fprintf(file,'\\begin{tabular}{@{}c|*{%d}{c}@{}}\n',length(kappas));
                    fprintf(file,...
                        mat2Table(corner,head,index',body,{'','%1.3e'},{},...
                        'headHook','\\toprule\\\\\n'));
                    fprintf(file,'\\end{tabular}\n\n\n');
                end
            end
            fprintf(file,'}\n');
        end
        if useEulerRef && ~(useExact && hasExact)
            latexCommand{end+1}=['\errDisEulerRef'];
            fprintf(file,'\\newcommand{%s}{\n',latexCommand{end});
            for tk=1:1:length(tEval)
                if mod(tEval(tk),0.1)==0
                    fprintf(file,...
                        'Error distributions at time %1.3f\n\n\n',tEval(tk));
                    corner='\\diagbox{Method}{$\\kappa$}';
                    head=kappas;
                    body=[];
                    index={};
                    if useEuler
                        temp1=errDisEulerRefEuler{tk};
                        temp=zeros(1,length(temp1));
                        for it1 = 1:1:length(temp1)
                            temp(it1)=mean(temp1{it1},1);
                        end
                        body=cat(1,body,temp);
                        index{end+1}='Euler';
                    end
                    if useMagnus1
                        temp1=errDisEulerRefMagnus1{tk};
                        temp=zeros(1,length(temp1));
                        for it1 = 1:1:length(temp1)
                            temp(it1)=mean(temp1{it1},1);
                        end
                        body=cat(1,body,temp);
                        index{end+1}='m1';
                    end
                    if useMagnus2
                        temp1=errDisEulerRefMagnus2{tk};
                        temp=zeros(1,length(temp1));
                        for it1 = 1:1:length(temp1)
                            temp(it1)=mean(temp1{it1},1);
                        end
                        body=cat(1,body,temp);
                        index{end+1}='m2';
                    end
                    if useMagnus3
                        temp1=errDisEulerRefMagnus3{tk};
                        temp=zeros(1,length(temp1));
                        for it1 = 1:1:length(temp1)
                            temp(it1)=mean(temp1{it1},1);
                        end
                        body=cat(1,body,temp);
                        index{end+1}='m3';
                    end
                    fprintf(file,'\\begin{tabular}{@{}c|*{%d}{c}@{}}\n',length(kappas));
                    fprintf(file,...
                        mat2Table(corner,head,index',body,{'','%1.3e'},{},...
                        'headHook','\\toprule\\\\\n'));
                    fprintf(file,'\\end{tabular}\n\n\n');
                end
            end
            fprintf(file,'}\n');
        end
        fclose(file);
    end
    function [latexFilePath,latexCommand]=saveComputationalTimes(saveAt)
        latexCommand={'\compTimes'};
        if strcmp(saveAt, '')
            latexFilePath='compTimes.tex';
        else
            latexFilePath=[saveAt,'/','compTimes.tex'];
            mkDir([tempPath,'/',saveAt]);
        end
        file = fopen([tempPath,'/',latexFilePath],'w+');
        if useExact
            fprintf(file,...
                    '\\newcommand{\\ctimeExact}{%g}\n',ctimeExactTotal);
        end
        if useEulerRef
            fprintf(file,...
                    '\\newcommand{\\ctimeEulerRef}{%g}\n',ctimeEulerRefTotal);
        end
        if useEuler
            fprintf(file,...
                    '\\newcommand{\\ctimeEuler}{%g}\n',ctimeEulerTotal);
        end
        if useMagnus1
        fprintf(file,...
                '\\newcommand{\\ctimeMagnusOne}{%g}\n',ctimeMagnus1Total);
        end
        if useMagnus2
        fprintf(file,...
                '\\newcommand{\\ctimeMagnusTwo}{%g}\n',ctimeMagnus2Total);
        end
        if useMagnus3
        fprintf(file,...
                '\\newcommand{\\ctimeMagnusThree}{%g}\n',ctimeMagnus3Total);
        end
        fprintf(file,...
                '\\newcommand{%s}{\n',latexCommand{1});
        fprintf(file,...
                '\\begin{tabular}{@{}*{3}{c}@{}}\n');
        fprintf(file,...
                'Method & Computational Time & Number of time evaluations\\\\\n');
        if useExact
            fprintf(file,...
                    'exact & $\\ctimeExact$ & $%d$\\\\\n',length(tnEuler)); 
        end
        if useEulerRef
            fprintf(file,...
                    'euler ref & $\\ctimeEulerRef$ & $%d$\\\\\n',length(tnEuler)); 
        end
        if useEuler
            fprintf(file,...
                    'euler & $\\ctimeEuler$ & $%d$\\\\\n',length(tnEuler)); 
        end
        if useMagnus1
            fprintf(file,...
                    'm1 & $\\ctimeMagnusOne$ & $%d$\\\\\n',length(tnMagnus)); 
        end
        if useMagnus2
            fprintf(file,...
                    'm2 & $\\ctimeMagnusTwo$ & $%d$\\\\\n',length(tnMagnus)); 
        end
        if useMagnus3
            fprintf(file,...
                    'm3 & $\\ctimeMagnusThree$ & $%d$\\\\\n',length(tnMagnus));
        end
        fprintf(file,...
                '\\end{tabular}\n');
        fprintf(file,...
                '}\n');    
        fclose(file);
    end
    function [latexFilePath,latexCommand]=saveModelCoefficients(saveAt)
        latexCommand={'\modelCoeff'};
        if strcmp(saveAt, '')
            latexFilePath='modelCoeff.tex';
        else
            latexFilePath=[saveAt,'/','modelCoeff.tex'];
            mkDir([tempPath,'/',saveAt]);
        end
        file = fopen([tempPath,'/',latexFilePath],'w+');
        fprintf(file,...
                '\\newcommand{\\funcH}{%s}\n',latex(sym(h)));
        fprintf(file,...
                '\\newcommand{\\funcFx}{%s}\n',latex(sym(fx)));
        fprintf(file,...
                '\\newcommand{\\funcFv}{%s}\n',latex(sym(fv)));
        fprintf(file,...
                '\\newcommand{\\funcGxx}{%s}\n',latex(sym(gxx)));
        fprintf(file,...
                '\\newcommand{\\funcGxv}{%s}\n',latex(sym(gxv)));
        fprintf(file,...
                '\\newcommand{\\funcGvv}{%s}\n',latex(sym(gvv)));
        fprintf(file,...
                '\\newcommand{\\funcSx}{%s}\n',latex(sym(sx)));
        fprintf(file,...
                '\\newcommand{\\funcSv}{%s}\n',latex(sym(sv)));
        fprintf(file,...
                '\\newcommand{\\funcPhi}{%s}\n',latex(sym(phi)));
        fprintf(file,...
                '\\newcommand{%s}{\n',latexCommand{1});
        fprintf(file,...
                '\\begin{align*}\n');
        fprintf(file,...
                'h(x,v)&\\coloneqq \\funcH\\\\\n');
        fprintf(file,...
                'f^x(x,v)&\\coloneqq \\funcFx\\\\\n');
        fprintf(file,...
                'f^v(x,v)&\\coloneqq \\funcFv\\\\\n');
        fprintf(file,...
                'g^{xx}(x,v)&\\coloneqq \\funcGxx\\\\\n');
        fprintf(file,...
                'g^{xv}(x,v)&\\coloneqq \\funcGxv\\\\\n');
        fprintf(file,...
                'g^{vv}(x,v)&\\coloneqq \\funcGvv\\\\\n');
        fprintf(file,...
                '\\sigma^{x}(x,v)&\\coloneqq \\funcSx\\\\\n');
        fprintf(file,...
                '\\sigma^{v}(x,v)&\\coloneqq \\funcSv\\\\\n');
        fprintf(file,...
                '\\phi(x,v)&\\coloneqq \\funcPhi\n');
        fprintf(file,...
                '\\end{align*}\n');
        fprintf(file,...
                '}\n');
        fclose(file);
    end
    function [latexFilePath,latexCommand]=saveNumericalParameters(saveAt)
        latexCommand={'\timeParam','\positionParam','\velocityParam'};
        if strcmp(saveAt, '')
            latexFilePath='numericalParam.tex';
        else
            latexFilePath=[saveAt,'/','numericalParam.tex'];
            mkDir([tempPath,'/',saveAt]);
        end
        file = fopen([tempPath,'/',latexFilePath],'w+');
        fprintf(file,...
                '\\newcommand{\\tUpper}{%g}\n',T);
        fprintf(file,...
                '\\newcommand{\\dtEulerRef}{%g}\n',dtEulerRef);
        fprintf(file,...
                '\\newcommand{\\dtEuler}{%g}\n',dtEuler);
        fprintf(file,...
                '\\newcommand{\\dtMagnus}{%g}\n',dtMagnus);
        fprintf(file,...
                '\\newcommand{\\xLower}{%g}\n',xl);
        fprintf(file,...
                '\\newcommand{\\xUpper}{%g}\n',xu);
        fprintf(file,...
                '\\newcommand{\\nx}{%g}\n',nx);
        fprintf(file,...
                '\\newcommand{\\vLower}{%g}\n',vl);
        fprintf(file,...
                '\\newcommand{\\vUpper}{%g}\n',vu);
        fprintf(file,...
                '\\newcommand{\\nv}{%g}\n',nv);
        fprintf(file,...
                '\\newcommand{\\simM}{%g}\n',M);
        fprintf(file,...
                '\\newcommand{%s}{\n',latexCommand{1});
        fprintf(file,...
                '\\begin{tabular}{@{}*{6}{c}@{}}\n');
        fprintf(file,...
                '$t_0$ & $T$ & time step $\\Delta^{\\text{euler ref}}$ & time step $\\Delta^{\\text{euler}}$ & time step $\\Delta^{\\text{magnus}}$ & simulations $M$ \\\\\n');
        fprintf(file,...
                '$0$ & $\\tUpper$ & $\\dtEulerRef$ & $\\dtEuler$ & $\\dtMagnus$ & $\\simM$ \\\\\n'); 
        fprintf(file,...
                '\\end{tabular}\n');
        fprintf(file,...
                '}\n');
        fprintf(file,...
                '\\newcommand{%s}{\n',latexCommand{2});
        fprintf(file,...
                '\\begin{tabular}{@{}*{3}{c}@{}}\n');
        fprintf(file,...
                'lower bound $a_x$ & upper bound $b_x$ & number of points $n_x$ \\\\\n');
        fprintf(file,...
                '$\\xLower$ & $\\xUpper$ & $\\nx$ \\\\\n'); 
        fprintf(file,...
                '\\end{tabular}\n');
        fprintf(file,...
                '}\n');
        fprintf(file,...
                '\\newcommand{%s}{\n',latexCommand{3});
        fprintf(file,...
                '\\begin{tabular}{@{}*{3}{c}@{}}\n');
        fprintf(file,...
                'lower bound $a_v$ & upper bound $b_v$ & number of points $n_v$ \\\\\n');
        fprintf(file,...
                '$\\vLower$ & $\\vUpper$ & $\\nv$ \\\\\n'); 
        fprintf(file,...
                '\\end{tabular}\n');
        fprintf(file,...
                '}\n');
        
        fclose(file);
    end
    function latexFilePath=saveFiguresLandscape(figures,saveAt,figName)
        if strcmp(figName,'')
            figName='fig';
        end
        if strcmp(saveAt, '')
            latexFilePath='figure.tex';
            latexFolderPath=tempPath;
        else
            latexFilePath=[saveAt,'/','figure.tex'];
            latexFolderPath=[tempPath,'/',saveAt];
            mkDir([tempPath,'/',saveAt]);
        end
        file = fopen([tempPath,'/',latexFilePath],'w+');
        for i=1:1:length(figures)
            picName=sprintf('%s_%d',figName,i);
            picPath = [latexFolderPath,'/',picName,'.',picType];
            if strcmp(saveAt, '')
                relPicPath=picName;
            else
                relPicPath=[saveAt,'/',picName];
            end
            exportgraphics(figures(i),picPath,'BackgroundColor',backgroundColor)
            fprintf(file,...
                '\\begin{landscape}\n');
            fprintf(file,...
                '\\includegraphics[width=.95\\columnwidth]{%s}\n',...
                changeSlash(relPicPath));
            fprintf(file,...
                '\\end{landscape}\n');
        end
        fclose(file);
    end
    function latexFilePath=saveFigures(figures,saveAt,figName)
        if strcmp(figName,'')
            figName='fig';
        end
        if strcmp(saveAt, '')
            latexFilePath='figure.tex';
            latexFolderPath=tempPath;
        else
            latexFilePath=[saveAt,'/','figure.tex'];
            latexFolderPath=[tempPath,'/',saveAt];
            mkDir([tempPath,'/',saveAt]);
        end
        file = fopen([tempPath,'/',latexFilePath],'w+');
        for i=1:1:length(figures)
            picName=sprintf('%s_%d',figName,i);
            picPath = [latexFolderPath,'/',picName,'.',picType];
            if strcmp(saveAt, '')
                relPicPath=picName;
            else
                relPicPath=[saveAt,'/',picName];
            end
            exportgraphics(figures(i),picPath,'BackgroundColor',backgroundColor)
            fprintf(file,...
                '\\begin{minipage}[c][][c]{\\linewidth}\n');
            fprintf(file,...
                '\\includegraphics[width=.95\\columnwidth]{%s}\n',...
                changeSlash(relPicPath));
            fprintf(file,...
                '\\end{minipage}\n');
        end
        fclose(file);
    end
    function latexFilePath=saveVideos(videos,saveAt)
        if strcmp(saveAt, '')
            latexFilePath='video.tex';
            latexFolderPath=tempPath;
        else
            latexFilePath=[saveAt,'/','video.tex'];
            latexFolderPath=[tempPath,'/',saveAt];
            mkDir([tempPath,'/',saveAt]);
        end
        file = fopen([tempPath,'/',latexFilePath],'w+');
        for i=1:1:length(videos)
            temp = strsplit(videos{i},'/');
            currVidName = [temp{end},'.pdf'];
            currVid = [latexFolderPath,'/',currVidName];
            copyfile([videos{i},'.pdf'],currVid);
            if strcmp(saveAt, '')
                relVidPath=currVidName;
            else
                relVidPath=[saveAt,'/',currVidName];
            end
            fprintf(file,...
                '\\begin{minipage}[c][][c]{\\linewidth}\n');
            fprintf(file,...
                '\\animategraphics[width=\\linewidth,autoplay,controls]{12}{%s}{}{}\n',...
                changeSlash(replace(relVidPath,'.pdf','')));
            fprintf(file,...
                '\\end{minipage}\n');
        end
        fclose(file);
    end

end
%% Auxiliary functions
function str=changeSlash(str)
    for i=1:1:length(str)
        if strcmp(str(i),'/')
            str(i)='/';
        end
    end
end
function delFile(file)
    if exist(file)
        delete(file);
    end
end
function delDir(dir)
    if exist(dir)==7
        rmdir(dir,'s');
    end
end
function cleanDir(mdir,except)
    except{end+1}='.';
    except{end+1}='..';
    for d = dir(mdir).'
      if ~any(strcmp(d.name,except))
          if d.isdir
              rmdir([d.folder,'/',d.name],'s');
          else
              delete([d.folder,'/',d.name]);
          end
      end
    end
end
function mkDir(dir)
    if exist(dir)==0
        mkdir(dir);
    end
end
function latexStr=mat2Table(corner,head,index,body,formatSpec,unit,varargin)
%     percent='\\,\\%%';
% formatSpec, unit order: index, body, head, corner
    headHook='';
    bodyHook='';
    for k=1:1:length(varargin)
        switch varargin{k}
            case 'headHook'
                headHook=varargin{k+1};
            case 'bodyHook'
                bodyHook=varargin{k+1};
        end
    end
    fSlen=length(formatSpec);
    tempFormatSpec={'','','',''};
    if fSlen>0 && fSlen < 4
        tempFormatSpec(1:fSlen)=formatSpec(:);
    end
    formatSpec=tempFormatSpec;
    ulen=length(unit);
    tempUnit={'','','',''};
    if ulen>0 && ulen < 4
        tempUnit(1:ulen)=unit(:);
    end
    unit=tempUnit;
    latexStr='';
    function x=functor(x,formatSpec,unit)
        if ~ischar(x)
            if mod(abs(x),1)==0 
                x=num2str(x,'%d');
            else
                x=num2str(x,formatSpec);
            end
        end
        x=[x,unit];
    end
    if ~isempty(corner)
        cornerStr=functor(corner,formatSpec{4},unit{4});
        latexStr=[latexStr,cornerStr,' & '];
    end
    if ~isempty(head)
        if iscell(head)
            headCell=cellfun(@(x)functor(x,formatSpec{3},unit{3}),head,...
                'UniformOutput',false);
        else
            headCell=arrayfun(@(x)functor(x,formatSpec{3},unit{3}),head,...
                'UniformOutput',false);
        end
        headStr=join(headCell,' & ');
        latexStr=[latexStr,headStr{1},'\\\\\n'];
    end
    latexStr=[latexStr,headHook];
    bodyCell=arrayfun(@(x)functor(x,formatSpec{2},unit{2}),body,...
        'UniformOutput',false);
    if ~isempty(index)
        if iscell(index)
            indexCell=cellfun(@(x)functor(x,formatSpec{1},unit{1}),index,...
                'UniformOutput',false);
        else
            indexCell=arrayfun(@(x)functor(x,formatSpec{1},unit{1}),index,...
                'UniformOutput',false);
        end
        bodyCell=cat(2,indexCell,bodyCell);
    end
    bodyStr=join(join(bodyCell,' & ',2),'\\\\\n',1);
    latexStr=[latexStr,bodyStr{1}];
    latexStr=[latexStr,bodyHook];
end