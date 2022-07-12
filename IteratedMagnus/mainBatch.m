% mainFunction(useEulerRef,useEuler,useMagnus1,useMagnus2,useMagnus3,useExact,methodName,...
%                       T,dT,dtEulerRef,dtEuler,nx,M)

% for nx=200, exact1
% dT=.05;
% for nx=300, exact1
% dT=.025;
% for nx=400, exact1
% dT=.01;
% for nx=200, Langevin1
% dT=.025;
% for nx=300, Langevin1
% dT=.01;
M=100;
for dtEulerRef=[1e-5]
    for dtEuler=[1e-4]
            for nx=[300]
				% const case
				% if nx<=100
					% dtMagnus = 1e-1;
				% elseif nx<=200
					% dtMagnus = 5e-2;
				% elseif nx<=300
					% dtMagnus = 2.5e-2;
				% else
					% dtMagnus = 1e-2;
				% end
				%  Langevin1
				if nx<=100
					dtMagnus = 5e-2;
				elseif nx<=200
					dtMagnus = 2.5e-2;
				elseif nx<=300
					dtMagnus = 1e-2;
				else
					dtMagnus = 5e-3;
				end
                fprintf('Start test with dtEulerRef=%1.3e, ',dtEulerRef)
                fprintf('dtEuler=%1.3e, ',dtEuler)
                fprintf('dtMagnus=%1.3e, ',dtMagnus)
                fprintf('M=%d, ',M)
                fprintf('nx=%d\n',nx)
                ticMainFunc = tic;
                mainFunction(1,1,1,1,1,1,'Langevin1',...
                             1,dtMagnus,dtEulerRef,dtEuler,nx,M)
                ctimeMainFunc = toc(ticMainFunc);
                fprintf('Elapsed time for all %1.3f\n',ctimeMainFunc)
                fprintf('###########################\n')
            end
    end
end