clear all; close all;
tempPath='ErrorLvl';
figPath=[pwd,'/',tempPath];
fileDir=[pwd,'/',tempPath,'/*.mat'];
files=dir(fileDir);

backgroundColor = [53,54,58]./255;
textColor = [237,237,237]./255;
backgroundColor = 'w';
textColor='k';

ref='Exact';
% errType='absError';
errType='relError';

saveDir=[figPath,'/',errType];
mkDir(saveDir);

excludeMethods={'Magnus1'};

errLvls=[1e-6,1e-5,1e-4,1e-3,1e-2,1e-1];%right bounds

ctimes=cell(length(errLvls),1);
dims=cell(length(errLvls),1);
dts=cell(length(errLvls),1);
errors=cell(length(errLvls),1);
methods=cell(length(errLvls),1);
methodsStr=cell(length(errLvls),1);
orders=cell(length(errLvls),1);
logDts=cell(length(errLvls),1);
for iF = 1:1:length(files)
    file=[files(iF).folder,'/',files(iF).name];
    curr=load(file);
    tempMethods=fieldnames(curr.Result);
    for iM = 1:1:length(tempMethods)
        method=tempMethods{iM};
        if any(contains(excludeMethods,method),'all')
            continue;
        end
        currError = curr.Result.(method).(ref).(errType)(end);
        currCtime = curr.Result.(method).ctime;
        currDt = curr.Result.(method).dt;
        currD = curr.Result.(method).nx;
        currMethod = replace(method,'Ref','');
        currMethodStr = [currMethod,sprintf(', dt=%1.1e, d=%d',currDt,currD)];
        currOrder = str2num(currMethod(end));
        if ~isempty(currOrder)
            currMethod=currMethod(1:end-1);
        else
            currOrder=-1;
        end
        currMethod = [currMethod(1),sprintf(', %1.1e, d=%d',currDt,currD)];
        currLogDt=-1;
        if strcmp(currMethod(1),'M')
            currLogDt = curr.Result.Euler.dt;
            currMethod = [currMethod,sprintf(', %1.1e',currLogDt)];
            currMethodStr = [currMethodStr,sprintf(', %1.1e',currLogDt)];
        end
        eLvl=find(currError<=errLvls,1,'first');
        if ~isempty(eLvl)
            ctimes{eLvl}(end+1)=currCtime;
            errors{eLvl}(end+1)=currError;
            methods{eLvl}{end+1}=currMethod;
            methodsStr{eLvl}{end+1}=currMethodStr;
            dims{eLvl}{end+1}=currD;
            dts{eLvl}{end+1}=currDt;
            orders{eLvl}{end+1}=currOrder;
            logDts{eLvl}{end+1}=currLogDt;
        end
    end
end

for iE = 1:1:length(errLvls)
    currMethods=methods{iE};
    currMethodsStr=methodsStr{iE};
    currCtimes=ctimes{iE};
    currErrors=errors{iE};
    currOrders=orders{iE};
    currDims=dims{iE};
    currDts=dts{iE};
    currLogDts=logDts{iE};
    if ~isempty(currMethods)
        fig=newFigure('visible','on',...
                      'backgroundColor',backgroundColor,...
                      'textColor',textColor);hold on;
        
        [uniqueMethodsStr,iU,~]=unique(currMethodsStr);
        x = currMethods(iU);
        c = reshape(currCtimes(iU),[],1);
        e = reshape(currErrors(iU),[],1);
        o = cell2mat(currOrders(iU));
        dt = cell2mat(currDts(iU));
        d = cell2mat(currDims(iU));
        logDt = cell2mat(currLogDts(iU));
        [x,c,e,o,xD,yD,xEuler,xMagnus2,xMagnus3,yEuler,yMagnus2,yMagnus3]=sortMethods(d,x,c,e,o,dt,logDt);
        [dUnique,ia,ic]=unique(d,'sorted');
        t = tiledlayout(1,length(x),'TileSpacing','none');
%         t = tiledlayout('flow','TileSpacing','none');
        bgAx = axes(t,'XTick',[],'YTick',[],'Box','off');
%         bgAx.Layout.TileSpan = [1 length(x)];
        ax = {};
%         currEX = categorical(x(xEuler),string(x(xEuler)));
%         currM3X = categorical(x(xMagnus3),string(x(xMagnus3)));
%         currM2X = categorical(x(xMagnus2),string(x(xMagnus2)));
%         currX=[currEX,currM3X,currM2X];
        oldTiles=1;
        for iD = 1:1:length(dUnique)
            hold off;
            xDim=xD==dUnique(iD);
            yDim=yD==dUnique(iD);
            currEX = categorical(x(xEuler & xDim),string(x(xEuler & xDim)));
            currM3X = categorical(x(xMagnus3 & xDim),string(x(xMagnus3 & xDim)));
            currM2X = categorical(x(xMagnus2 & xDim),string(x(xMagnus2 & xDim)));
            currX=[currEX,currM3X,currM2X];
            currTiles=sum(xEuler & xDim,'all')+sum((xMagnus3 & xDim) | (xMagnus2 & xDim),'all');
%             ax{end+1} = axes(t);
            ax{end+1} = nexttile(oldTiles,[1 currTiles]);hold on;
            set(gca, 'color', backgroundColor);
            set(gca, 'XColor', textColor);
%             set(gca,'FontSize',22)
            currAx=ax{end};
%             currAx.Layout.Tile = iD;
%             currAx.Layout.TileSpan = [1 currTiles];hold on;
            oldTiles=oldTiles+currTiles;
            currAx.Box = 'off';
            title(currAx,sprintf('Spatial dim: %d',dUnique(iD)),...
                'Interpreter','latex',...
                'Color',textColor,...
                'FontSize',16);
            
            yyaxis left
            % ctimes euler
            if length(currEX)>0
                bEuler=ctimeBars(currEX,c(yEuler & yDim),[],textColor,[0 0.4470 0.7410],'bottom');
            end
            % ctimes Magnus 3
            if length(currM3X)>0
                bMagnus3=ctimeBars(currM3X,c(yMagnus3 & yDim),[],textColor,[0 0.4470 0.7410],'bottom');
            end
            % ctimes Magnus 2
            if length(currM2X)>0
                bMagnus2=ctimeBars(currM2X,c(yMagnus2 & yDim),.25,textColor,[0.3010 0.7450 0.9330],'top');
            end
            set(gca, 'YScale', 'log')
            tempYlim=ylim;
            tempYlim=tempYlim + tempYlim/3 .* [-1 1];
            ylim(tempYlim)
            if iD>1
                currAx.YAxis(1).TickValues = [];
                currAx.YAxis(1).TickLabels = {};
                currAx.YAxis(1).LineWidth = 1;
                currAx.YAxis(1).Color = textColor;
%                 currAx.YAxis(1).Visible = 'off';
            end
            if iD==1
                ylabel('Computational time in s')
            end
            yyaxis right
            % errors euler
            if length(currEX)>0
                bEuler=errorBars(currEX,e(yEuler & yDim),[],textColor,[0.6350 0.0780 0.1840],'bottom');
            end
            % errors Magnus 3
            if length(currM3X)>0
                bMagnus3=errorBars(currM3X,e(yMagnus3 & yDim),[],textColor,[0.6350 0.0780 0.1840],'bottom');
            end
            % errors Magnus 2
            if length(currM2X)>0
                bMagnus2=errorBars(currM2X,e(yMagnus2 & yDim),.25,textColor,[0.8500 0.3250 0.0980],'top');
            end
            tempYlim=ylim;
            tempYlim=tempYlim + tempYlim/5 .* [-1 1];
            ylim(tempYlim)
            if iD<length(dUnique)
%                 currAx.YAxis(2).Visible = 'off';
                currAx.YAxis(2).TickValues = [];
                currAx.YAxis(2).TickLabels = {};
                currAx.YAxis(2).LineWidth = 1;
                currAx.YAxis(2).Color = textColor;
            else
                ylabel('Error')
            end
            XTickLabelRotation=90;
            currXTickLabels=xticklabels;
            xticklabels(replace(replace(currXTickLabels,[', d='+digitsPattern],' '),'*',''));
%             currAx.XAxis.TickLabelInterpreter='latex';
%             currAx.XAxis.TickLabelColor=textColor;
%             xlabel(replace(replace(currXTickLabels,[', d='+digitsPattern],' '),'*','\\newline'),"Interpreter","latex")
            currAx.XAxis.Color=textColor;
            currAx.XAxis.LineWidth=3;
        end
        leftYAxis={};
        rightYAxis={};
        leftMax=[1e6,0];
        rightMax=[1e6,0];
        for iA=1:1:length(ax)
            cax=ax{iA};
            leftYAxis{end+1}=cax.YAxis(1);
            rightYAxis{end+1}=cax.YAxis(2);
            leftMax(1)=min([leftMax(1),cax.YAxis(1).Limits(1)]);
            leftMax(2)=max([leftMax(2),cax.YAxis(1).Limits(2)]);
            rightMax(1)=min([rightMax(1),cax.YAxis(2).Limits(1)]);
            rightMax(2)=max([rightMax(2),cax.YAxis(2).Limits(2)]);
        end
        for iA=1:1:length(ax)
            cax=ax{iA};
            cax.YAxis(1).Limits=leftMax;
            cax.YAxis(2).Limits=rightMax;
        end
        leftYAxis=[leftYAxis{:}];
        rightYAxis=[rightYAxis{:}];
        linkprop(leftYAxis(:),'Limits');
        linkprop(rightYAxis(:),'Limits');
        title(t,sprintf('Error level: $\\left(%1.1e,%1.1e\\right]$',...
                errLvls(iE)/10,errLvls(iE)),...
                'Interpreter','latex',...
                'Color',textColor,...
                'FontSize',22);
        exportgraphics(fig,[saveDir,'/',sprintf('errorlvl_%1.1e.pdf',errLvls(iE))],...
                       "BackgroundColor",backgroundColor);
    end
end
function b=ctimeBars(x,c,w,textColor,faceColor,vAlignment)
if isempty(w)
    b=bar(x,[c,nan.*ones(length(x),1)],'FaceColor',faceColor);
else
    b=bar(x,[c,nan.*ones(length(x),1)],w,'FaceColor',faceColor);
end
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = strsplit(sprintf("%1.1f ",b(1).YData),' ');
text(xtips1,ytips1,labels1(1:end-1),'HorizontalAlignment','center',...
    'VerticalAlignment',vAlignment,'Color',textColor)
end
function b=errorBars(x,e,w,textColor,faceColor,vAlignment)
if isempty(w)
    b=bar(x,[nan.*ones(length(x),1),e],'FaceColor',faceColor);
else
    b=bar(x,[nan.*ones(length(x),1),e],w,'FaceColor',faceColor);
end
xtips2 = b(2).XEndPoints;
ytips2 = b(2).YEndPoints;
labels2 = strsplit(sprintf("%1.2e ",b(2).YData),' ');
text(xtips2,ytips2,labels2(1:end-1),'HorizontalAlignment','center',...
    'VerticalAlignment',vAlignment,'Color',textColor)
end
function [x,c,e,o,xD,yD,xEuler,xMagnus2,xMagnus3,yEuler,yMagnus2,yMagnus3]=...
    sortMethods(d,x,c,e,o,dt,logDt)
    [d,iSorted] = sort(d);
    x=x(iSorted);
%     x=categorical(x(iSorted),'Ordinal',true);
    c=c(iSorted);
    e=e(iSorted);
    o=o(iSorted);
    dt=dt(iSorted);
    logDt=logDt(iSorted);
    
    paramArray = [o',dt',d',logDt'];

    [~,ia,~] = unique(paramArray(:,2:end),'rows','stable');
    
    otemp=paramArray(ia,1)';
    x=x(ia);
    xD=d(ia);
    yD=d;
    xEuler=-1==otemp;
    xMagnus2=2==otemp;
    xMagnus3=(3==otemp | 2==otemp); % could make trouble if orders are not ordered in paramArray, if necessary first order them
    yEuler=-1==o;
    yMagnus2=2==o;
    yMagnus3=3==o;

    x=categorical(x,x,'Ordinal',true);
end
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