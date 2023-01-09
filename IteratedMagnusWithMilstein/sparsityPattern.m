function patterns=sparsityPattern(A,B,numPoints,varargin)
showTitle=true;
backgroundColor='w';
textColor='k';
for iV=1:2:length(varargin)
    switch varargin{iV}
        case 'showTitle'
            showTitle=varargin{iV+1};
        case 'backgroundColor'
            backgroundColor = varargin{iV+1};
        case 'textColor'
            textColor = varargin{iV+1};
    end
end
% nnz(A)/prod(size(A))
% nnz(B)/prod(size(B))
% whos('A')
% numPoints=100;
d=size(A,1);
c=floor(d/2);
lb=1;
ub=d;
xStr='';
if d>numPoints
    lb=c-floor(numPoints/2);
    ub=c+ceil(numPoints/2);
    xStr=sprintf('Central %d x %d points of entire matrix',min(d,numPoints),min(d,numPoints));
end
patterns={};
[plots,legendEntries]=beginFigure();
yStr=sprintf('In total %d non-zero diagonals',size(spdiags(A),2));
spy(A(lb:ub,lb:ub));xticklabels('');yticklabels('');
endFigure(plots,legendEntries,xStr,yStr,'Sparsity pattern of $A$')

[plots,legendEntries]=beginFigure();
yStr=sprintf('In total %d non-zero diagonals',size(spdiags(B),2));
spy(B(lb:ub,lb:ub));xticklabels('');yticklabels('');
endFigure(plots,legendEntries,xStr,yStr,'Sparsity pattern of $B$')

BA=B*A-A*B;
[plots,legendEntries]=beginFigure();
yStr=sprintf('In total %d non-zero diagonals',size(spdiags(BA),2));
spy(BA(lb:ub,lb:ub));xticklabels('');yticklabels('');
endFigure(plots,legendEntries,xStr,yStr,'Sparsity pattern of $\left[B,A\right]$')

BAA=BA*A-A*BA;
[plots,legendEntries]=beginFigure();
yStr=sprintf('In total %d non-zero diagonals',size(spdiags(BAA),2));
spy(BAA(lb:ub,lb:ub));xticklabels('');yticklabels('');
endFigure(plots,legendEntries,xStr,yStr,'Sparsity pattern of $\left[\left[B,A\right],A\right]$')

BAB=BA*B-B*BA;
[plots,legendEntries]=beginFigure();
yStr=sprintf('In total %d non-zero diagonals',size(spdiags(BAB),2));
spy(BAB(lb:ub,lb:ub));xticklabels('');yticklabels('');
endFigure(plots,legendEntries,xStr,yStr,'Sparsity pattern of $\left[\left[B,A\right],B\right]$')
patterns=cat(2,patterns{:});
function [plots,legendEntries]=beginFigure()
    patterns{end+1}=newFigure('backgroundColor',backgroundColor,'textColor',textColor);
    legendEntries = {};
    plots = [];
end
function endFigure(plots,legendEntries,xLabel,yLabel,titleStr)
    if ~isempty(legendEntries)
        legend(plots,legendEntries,...
          'Location','southoutside',...
          'NumColumns',2,...
          'Interpreter','latex',...
          'TextColor',textColor); 
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
