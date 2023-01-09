function figures=plotErrorDis(errDisRefApprox,xRefApprox,...
                              fRefApprox,nBlowupsRefApprox,...
                              kappas,tEval,varargin)
%%PLOTERRORDIS plots error distribution. 
%   Input:
%   Output:
%   Usage:
%
% See also exact, magnusConst, eulerConst.
showTitle=true;
normType='rel1';
bw = 1;
percentile=1;
fixScale=true;
backgroundColor='w';
textColor='k';
for iV=1:2:length(varargin)
    switch varargin{iV}
        case 'showTitle'
            showTitle=varargin{iV+1};
        case 'normType'
            normType=varargin{iV+1};
        case 'blackWhite'
            bw=double(1+varargin{iV+1});
        case 'percentile'
            percentile=varargin{iV+1};
        case 'fixScale'
            fixScale=varargin{iV+1};
        case 'backgroundColor'
            backgroundColor = varargin{iV+1};
        case 'textColor'
            textColor = varargin{iV+1};
    end
end


figures=cell(length(tEval),length(errDisRefApprox{1}));

[exact,eulerRef,euler,m1,m2,m3,m4,inputNameDict,colorDict]=methodColors();
methods=strsplit(regexprep(replace(inputname(1),'errDis',''), '([A-Z])', ' $1'),' ');
i=2;
tempName= methods{2};
if strcmp(methods{2},'Euler')
if strcmp(methods{3},'Ref')
    tempName=[tempName,methods{3}];
    i=i+1;
end
end
refName=inputNameDict.(tempName);
i=i+1;
tempName= methods{i};
if strcmp(methods{i},'Euler')
if (i<length(methods)) && (strcmp(methods{i+1},'Ref'))
    tempName=[tempName,methods{i+1}];
end
end
approxName=inputNameDict.(tempName);
approxColor=colorDict.(tempName);
% refName=methodDict.(methods{2});
% approxName=methodDict.(methods{3});
%% Figure 1
for iKappa=1:1:length(kappas)
    if fixScale
        temp=[xRefApprox{:}];
        xmin=1e6;
        xmax=-1;
        for tj=1:1:size(temp,2)
            xmin=min(xmin,min(temp{iKappa,tj},[],'all'));
            xmax=max(xmax,max(temp{iKappa,tj},[],'all'));
        end
%         temp=cell2mat(temp(iKappa,:));
%         xmin=min(temp,[],'all');
%         xmax=max(temp,[],'all');
    end
    for tk=1:1:length(tEval)
        fig=newFigure('backgroundColor',backgroundColor,'textColor',textColor);
        [plots,legendEntries]=beginFigure(tk,iKappa);
        if fixScale
            xlim([xmin,xmax]);
        end
        plot_CDF(tk,iKappa)
        endFigure(plots,legendEntries,'x',sprintf('CDF_{Err_{%2.3g}}(x)',tEval(tk)),...
            sprintf('Empirical CDF of %s vs %s at t=%1.3g for kappa=%d',...
            refName,approxName,tEval(tk),kappas(iKappa)));
    end
end

%% Plot functions
function plot_CDF(tk,iKappa)
    plots(end+1)=plot(xRefApprox{tk}{iKappa},fRefApprox{tk}{iKappa},...
                    'LineStyle',...
                    approxColor.linestyle{bw},...
                    'Marker',...
                    approxColor.marker{bw},...
                    'Color',...
                    approxColor.color{bw});
       if nBlowupsRefApprox{tk}{iKappa}>-1
           xLim=xlim;
           yLim=ylim;
           text(xLim(2),yLim(1),sprintf('Number of excluded blowups %d',nBlowupsRefApprox{tk}{iKappa}),...
               'HorizontalAlignment','right','VerticalAlignment','bottom','Fontsize',22,'Color',textColor);
       end
       text(xLim(1),yLim(2),sprintf('E[Err_t]=%1.3e',mean(errDisRefApprox{tk}{iKappa},'all')),...
               'HorizontalAlignment','left','VerticalAlignment','top','Fontsize',22,'Color',textColor);
       legendEntries{end+1}=approxName;
end

%% Auxiliary functions
function [plots,legendEntries]=beginFigure(tk,iKappa)
    figures{tk,iKappa}=fig;
    legendEntries = {};
    plots = [];
end
function endFigure(plots,legendEntries,xLabel,yLabel,titleStr)
    legend(plots,legendEntries,...
      'Location','southoutside',...
      'NumColumns',4,...
      'Interpreter','latex',...
      'TextColor',textColor); 
    switch normType
        case {'rel1','rel2'}
            ytickformat('%1.2f')
            a=[cellstr(num2str(get(gca,'xtick')'*100))]; 
            pct = char(ones(size(a,1),1)*'%');
            new_xticks = [char(a),pct];
            set(gca,'xticklabel',new_xticks)
            a=[cellstr(num2str(get(gca,'ytick')'))]; 
            a{end}=char(num2str(percentile));
            new_yticks = [char(a)];
            set(gca,'yticklabel',new_yticks)
            xlabel(xLabel, 'fontweight', 'bold')
            ylabel(yLabel, 'fontweight', 'bold')
            title(titleStr,'Interpreter','latex')
%         case 'abs1'
    end
    if ~strcmp(xLabel,'')
        xlabel(xLabel, 'fontweight', 'bold')
    end
    if ~strcmp(yLabel,'')
        ylabel(yLabel, 'fontweight', 'bold')
    end
    if ~strcmp(titleStr,'') && showTitle
        title(titleStr,'Color',textColor,'Interpreter','latex')
    end
end
end