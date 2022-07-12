clear all; close all;
tempPath='LogTest';
figPath=[pwd,'/',tempPath];
fileDir=[pwd,'/',tempPath,'/*.mat'];
files=dir(fileDir);

% backgroundColor = [53,54,58]./255;
% textColor = [237,237,237]./255;
backgroundColor = 'None';
textColor='k';

ref='Exact';
% errType='absError';
errType='relError';

saveDir=[figPath,'/',errType];
mkDir(saveDir);

excludeMethods={};


ctimes=zeros(length(files),3);
errors=zeros(length(files),3);
logDt=zeros(length(files),3);

for iF = 1:1:length(files)
    file=[files(iF).folder,'/',files(iF).name];
    curr=load(file);
    tempMethods=fieldnames(curr.logTest);
    for iM = 1:1:length(tempMethods)
        method=tempMethods{iM};
        if any(contains(excludeMethods,method),'all')
            continue;
        end   
        ctimes(iF,str2num(method(end)))=curr.logTest.(method).ctime;
        errors(iF,str2num(method(end)))=curr.logTest.(method).(ref).(errType)(end);
        logDt(iF,str2num(method(end)))= curr.logTest.(method).dtLog;
    end
end
M=curr.logTest.Magnus3.M;
dt=curr.logTest.Magnus3.dt;
d=curr.logTest.Magnus3.nx;
T=curr.logTest.Magnus3.T;

[logDtSorted,I]=sort(logDt,1);
ctimesSorted=ctimes(I(:,1),:)./M;
if strcmp(errType,'relError')
    errorsSorted=errors(I(:,1),:).*100;
else
    errorsSorted=errors(I(:,1),:);
end
errorsSorted(isnan(errorsSorted))=1e2;

legendEntry={};
fig=newFigure('visible','on',...
              'backgroundColor',backgroundColor,...
              'textColor',textColor);hold on;

title(sprintf('T=%1.2f, M=%d, d=%d, dt=%1.2e',T,M,d,dt),...
        'Interpreter','latex',...
        'Color',textColor,...
        'FontSize',16);
xlabel('$\Delta_t^{\mathrm{Leb}}$','Interpreter','latex')
set(gca, 'XScale', 'log')

yyaxis left; 
plot(logDtSorted(:,2),ctimesSorted(:,2),'Color',[0.3010 0.7450 0.9330],'LineStyle',':')
legendEntry{end+1}='ctime m2';
plot(logDtSorted(:,3),ctimesSorted(:,3),'Color',[0 0.4470 0.7410],'LineStyle','--')
legendEntry{end+1}='ctime m3';
ylabel('Avg. Comp. Time in sec. per simulation (log scale)')
set(gca, 'YScale', 'log')

yyaxis right; 
plot(logDtSorted(:,2),errorsSorted(:,2),'Color',[0.8500 0.3250 0.0980],'LineStyle',':')
legendEntry{end+1}='error m2';
plot(logDtSorted(:,3),errorsSorted(:,3),'Color',[0.6350 0.0780 0.1840],'LineStyle','--')
legendEntry{end+1}='error m3';
if strcmp(errType,'relError')
    ylabel('Error in % (log scale)')
else
    ylabel('Error (log scale)')
end
set(gca, 'YScale', 'log')
p=prctile(errors,.99,1);
yLimR=ylim;
yLimR(2)=max(p(2:3))+.5*yLimR(1);
ylim(yLimR)

legend(legendEntry,'Location','southoutside','NumColumns',4,'Color',backgroundColor)
exportgraphics(fig,[saveDir,'/',sprintf('LogTest_T%1.2f_M%d_d%d_dt%1.2e.pdf',T,M,d,dt)],...
                       "BackgroundColor",backgroundColor);
% dark blue: [0 0.4470 0.7410]
% orange: [0.8500 0.3250 0.0980]
% yellow: [0.9290 0.6940 0.1250]
% purple: [0.4940 0.1840 0.5560]
% green: [0.4660 0.6740 0.1880]
% light blue: [0.3010 0.7450 0.9330]
% red: [0.6350 0.0780 0.1840]
function delDir(dir)
    if exist(dir)==7
        rmdir(dir,'s');
    end
end
function mkDir(dir)
    if exist(dir)==0
        mkdir(dir);
    end
end
function delFile(file)
    if exist(file)
        delete(file);
    end
end