function figures=plotTrajectories(tEulerRef,utEulerRef,tEuler,utEuler,...
                          tMagnus,utMagnus1,utMagnus2,utMagnus3,d1,d2)
%%PLOTTRAJECTORIES
%
%

% exact.linestyle={'-','-'};
% exact.marker={'none','none'};
% exact.color={'k','k'};

if length(tEulerRef)<=length(tMagnus)
eulerRef.linestyle={'none','none'};
eulerRef.marker={'+','+'};  
else
eulerRef.linestyle={'-','-'};
eulerRef.marker={'none','none'};
end
eulerRef.color={'k',rgb2gray([1 0 0])};

if length(tEuler)<=length(tMagnus)
euler.linestyle={'none','none'};
euler.marker={'*','*'};  
else
euler.linestyle={'--','--'};
euler.marker={'none','none'};
end
euler.color={'r',rgb2gray([1 0 0])};

m1.linestyle={'none','none'};
m1.marker={'.','.'};
m1.color={'b',rgb2gray([0 0 1])};

m2.linestyle={'none','none'};
m2.marker={'o','o'};
m2.color={'y',rgb2gray([1 1 0])};

m3.linestyle={'none','none'};
m3.marker={'x','x'};
m3.color={'g',rgb2gray([0 1 0])};

bw = 1;
figures=[];

% [~,w]=max(sum(abs(diff(utEuler(d1,d2,:,:),1,3)),3),[],4);
% fprintf('We used trajectory %d for the plot.\n',w)
w=1;

[plots,legendEntries]=beginFigure();
plotEulerRef();
plotEuler();
plotMagnus3();
endFigure(plots,legendEntries,...
          'Time',...
          sprintf('$u_t(x_i,v_j)$ at $i=%d$, $j=%d$',d1,d2))
[plots,legendEntries]=beginFigure();
plotEulerRef();
plotEuler();
plotMagnus2();
plotMagnus3();
endFigure(plots,legendEntries,...
          'Time',...
          sprintf('$u_t(x_i,v_j)$ at $i=%d$, $j=%d$',d1,d2))
[plots,legendEntries]=beginFigure();
plotEulerRef();
plotEuler();
plotMagnus1();
plotMagnus2();
plotMagnus3();
endFigure(plots,legendEntries,...
          'Time',...
          sprintf('$u_t(x_i,v_j)$ at $i=%d$, $j=%d$',d1,d2))
    

function [plots,legendEntries]=beginFigure()
    fig=figure('units','normalized',...
              'outerposition',[0 0 1 1]); hold on;
    fig.WindowState = 'minimized';
    figure_properties(fig);
    figures(end+1)=fig;
    legendEntries = {};
    plots = [];
end
function endFigure(plots,legendEntries,xLabel,yLabel)
    legend(plots,legendEntries,...
      'Location','southoutside',...
      'NumColumns',4,...
      'Interpreter','latex'); 
    xlabel(xLabel)
    ylabel(yLabel,'Interpreter','latex')
end
function plotEulerRef()
    plots(end+1)=plot(tEulerRef,squeeze(utEulerRef(d1,d2,:,w)),...
                    'LineStyle',...
                    eulerRef.linestyle{bw},...
                    'Marker',...
                    eulerRef.marker{bw},...
                    'Color',...
                    eulerRef.color{bw});
       legendEntries{end+1}='euler ref';
end
function plotEuler()
    plots(end+1)=plot(tEuler,squeeze(utEuler(d1,d2,:,w)),...
                    'LineStyle',...
                    euler.linestyle{bw},...
                    'Marker',...
                    euler.marker{bw},...
                    'Color',...
                    euler.color{bw});
       legendEntries{end+1}='euler';
end
function plotMagnus1()
    if ~isempty(utMagnus1)
       plots(end+1)=plot(tMagnus,squeeze(utMagnus1(d1,d2,:,w)),...
                    'LineStyle',...
                    m1.linestyle{bw},...
                    'Marker',...
                    m1.marker{bw},...
                    'Color',...
                    m1.color{bw});
       legendEntries{end+1}='m1';
    end
end
function plotMagnus2()
    if ~isempty(utMagnus1)
       plots(end+1)=plot(tMagnus,squeeze(utMagnus2(d1,d2,:,w)),...
                    'LineStyle',...
                    m2.linestyle{bw},...
                    'Marker',...
                    m2.marker{bw},...
                    'Color',...
                    m2.color{bw});
       legendEntries{end+1}='m2';
    end
end
function plotMagnus3()
    if ~isempty(utMagnus1)
       plots(end+1)=plot(tMagnus,squeeze(utMagnus3(d1,d2,:,w)),...
                    'LineStyle',...
                    m3.linestyle{bw},...
                    'Marker',...
                    m3.marker{bw},...
                    'Color',...
                    m3.color{bw});
       legendEntries{end+1}='m3';
    end
end
end
function figure_properties(fig)
fontsize=22;
linewidth=2;
markersize=12;
set(gca,'FontSize',fontsize)
set(gca,'defaultLineMarkerSize',markersize)
set(fig,'defaultlinelinewidth',linewidth)
set(fig,'defaultaxeslinewidth',linewidth)
set(fig,'defaultpatchlinewidth',linewidth)
set(fig,'defaultAxesFontSize',fontsize)
end