function figures=plotCDF(x,f,approxMethod,T,percentile,titleStr,varargin)
blowup=-1;
if ~isempty(varargin)
    for k=1:2:length(varargin)
        switch varargin{k}
            case 'blowup'
                blowup=varargin{k+1};
        end
    end
end
[exact,eulerRef,euler,m1,m2,m3,m4,inputNameDict]=methodColors();

bw = 1;
figures=[];

switch approxMethod
    case 'EulerRef'
        [plots,legendEntries]=beginFigure();
        plotEulerRef();
        endFigure(plots,legendEntries,...
                  'x',...
                  sprintf('CDF_{Err_{%2.3g}}(x)',T))
    case 'Euler'
        [plots,legendEntries]=beginFigure();
        plotEuler();
        endFigure(plots,legendEntries,...
                  'x',...
                  sprintf('CDF_{Err_{%2.3g}}(x)',T))
    case 'Magnus1'
        [plots,legendEntries]=beginFigure();
        plotMagnus1();
        endFigure(plots,legendEntries,...
                  'x',...
                  sprintf('CDF_{Err_{%2.3g}}(x)',T))
    case 'Magnus2'
        [plots,legendEntries]=beginFigure();
        plotMagnus2();
        endFigure(plots,legendEntries,...
                  'x',...
                  sprintf('CDF_{Err_{%2.3g}}(x)',T))
    case 'Magnus3'
        [plots,legendEntries]=beginFigure();
        plotMagnus3();
        endFigure(plots,legendEntries,...
                  'x',...
                  sprintf('CDF_{Err_{%2.3g}}(x)',T))
    
end
function [plots,legendEntries]=beginFigure()
    fig=newFigure();
    figures(end+1)=fig;
    legendEntries = {};
    plots = [];
end
function endFigure(plots,legendEntries,xLabel,yLabel)
    legend(plots,legendEntries,...
      'Location','southoutside',...
      'NumColumns',4,...
      'Interpreter','latex'); 
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
end
function plotEulerRef()
    if ~isempty(x)
       plots(end+1)=plot(x,f,...
                    'LineStyle',...
                    eulerRef.linestyle{bw},...
                    'Marker',...
                    eulerRef.marker{bw},...
                    'Color',...
                    eulerRef.color{bw});
       if blowup>-1
           xLim=xlim;
           yLim=ylim;
           text(xLim(2),yLim(1),sprintf('Number of excluded blowups %d',blowup),...
               'HorizontalAlignment','right','VerticalAlignment','bottom','Fontsize',22);
       end
       legendEntries{end+1}=inputNameDict.EulerRef;
    end
end
function plotEuler()
    if ~isempty(x)
       plots(end+1)=plot(x,f,...
                    'LineStyle',...
                    euler.linestyle{bw},...
                    'Marker',...
                    euler.marker{bw},...
                    'Color',...
                    euler.color{bw});
       if blowup>-1
           xLim=xlim;
           yLim=ylim;
           text(xLim(2),yLim(1),sprintf('Number of excluded blowups %d',blowup),...
               'HorizontalAlignment','right','VerticalAlignment','bottom','Fontsize',22);
       end
       legendEntries{end+1}=inputNameDict.Euler;
    end
end
function plotMagnus1()
    if ~isempty(x)
       plots(end+1)=plot(x,f,...
                    'LineStyle',...
                    m1.linestyle{bw},...
                    'Marker',...
                    m1.marker{bw},...
                    'Color',...
                    m1.color{bw});
       if blowup>-1
           xLim=xlim;
           yLim=ylim;
           text(xLim(2),yLim(1),sprintf('Number of excluded blowups %d',blowup),...
               'HorizontalAlignment','right','VerticalAlignment','bottom','Fontsize',22);
       end
       legendEntries{end+1}=inputNameDict.Magnus1;
    end
end
function plotMagnus2()
    if ~isempty(x)
       plots(end+1)=plot(x,f,...
                    'LineStyle',...
                    m2.linestyle{bw},...
                    'Marker',...
                    m2.marker{bw},...
                    'Color',...
                    m2.color{bw});
       if blowup>-1
           xLim=xlim;
           yLim=ylim;
           text(xLim(2),yLim(1),sprintf('Number of excluded blowups %d',blowup),...
               'HorizontalAlignment','right','VerticalAlignment','bottom','Fontsize',22);
       end
       legendEntries{end+1}=inputNameDict.Magnus2;
    end
end
function plotMagnus3()
    if ~isempty(x)
       plots(end+1)=plot(x,f,...
                    'LineStyle',...
                    m3.linestyle{bw},...
                    'Marker',...
                    m3.marker{bw},...
                    'Color',...
                    m3.color{bw});
       if blowup>-1
           xLim=xlim;
           yLim=ylim;
           text(xLim(2),yLim(1),...
               sprintf('Number of excluded blowups %d',blowup),...
               'HorizontalAlignment','right','VerticalAlignment','bottom','Fontsize',22);
       end
       legendEntries{end+1}=inputNameDict.Magnus3;
    end
end
end