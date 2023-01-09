function fileNames=videoErrorDis(gifPath,gifName,errorDisRefFig,varargin)
backgroundColor='w';
for iV=1:2:length(varargin)
    switch varargin{iV}
        case 'backgroundColor'
            backgroundColor = varargin{iV+1};
    end
end
mkDir(gifPath);
fileName=[gifPath,'/',gifName];
delete([fileName,'_*.*']);
[nT,nKappa]=size(errorDisRefFig);
for iKappa = 1:1:nKappa
    fileNames{iKappa}=[fileName,sprintf('_%d',iKappa)];
    for tk = 1:1:nT
        exportgraphics(errorDisRefFig{tk,iKappa},...
                       [fileName,sprintf('_%d',iKappa),'.gif'],...
                       'Append',true,...
                       'BackgroundColor', backgroundColor);
        exportgraphics(errorDisRefFig{tk,iKappa},...
                      [fileName,sprintf('_%d',iKappa),'.pdf'],...
                       'Append',true,...
                       'BackgroundColor', backgroundColor);
    end
end
end
function mkDir(dir)
    if exist(dir)==0
        mkdir(dir);
    end
end