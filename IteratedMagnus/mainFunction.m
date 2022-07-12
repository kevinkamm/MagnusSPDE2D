function mainFunction(useEulerRef,useEuler,useMagnus1,useMagnus2,useMagnus3,useExact,method,...
                      T,dT,dtEulerRef,dtEuler,nx,M)
% clear all; 
close all; fclose('all'); 
rng(0);
availGPU=gpuDeviceCount("available");
if availGPU > 0
    device='gpu';
    gpu=gpuDevice(1);
else
    device='cpu';
end
% don't change, need multiprocessing for coefficients and multithreading
% for Magnus
delete(gcp('nocreate'));parpool('threads'); % multiprocessing
%% Select methods
% choose methods to be computed
% useEulerRef = true;
% useEuler = true;
% useMagnus1 = true;
% useMagnus2 = true;
% useMagnus3 = true;
% useExact = true; %only use for constant coefficient case with set initial datum
useExactApprox = false;%only use for constant coefficient case
% select reference method for errors
% referenceMethod = 'exact';
% referenceMethod = 'eulerRef';
%% Error regions
kappas=[1,2,4];
blowup=1;
percentile=.99;
errorDisType='rel1';
%% Parameters
% Model
% constant coefficient case with exact solution
a=1.1/10*10;
sigma=1/10*sqrt(10);
[h,fx,fv,gxx,gxv,gvv,sx,sv,phi,hasExact,methodName]=...
    selectCoefficients(method,[a,sigma]);
% [h,fx,fv,gxx,gxv,gvv,sx,sv,phi,hasExact,methodName]=...
%     selectCoefficients('Langevin1',[a,sigma]);
%% Video & Plot parameters
makeErrorGif=false;
makePlots=false;
dirAbsErrorGif=['Results/Videos/',methodName,'/','AbsErrors'];
dirErrorDisGif=['Results/Videos/',methodName,'/','ErrorDis'];
% backgroundColor = 'k';
backgroundColor = [53,54,58]./255;
% textColor = 'w';
textColor = [237,237,237]./255;
%%
% Time
% T=1;
% T=.5;
% T=.1;
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

% dtEulerRef=1e-4;
% dtEuler=1e-3;
dtMagnus=dtEuler;
NeulerRef = ceil(dT/dtEulerRef)+1;
Neuler = ceil(dT/dtEuler)+1;
Nmagnus= ceil(dT/dtMagnus)+1;

% Simulations
% M = 100;
%% Space grids
% Position
xl=-4;
xu=4;
% nx=300;
xgrid=linspace(xl,xu,nx+2)'; %nx+2 x 1
% Velocity
nv=nx;
vl=xl;
vu=xu;
vgrid=linspace(vl,vu,nv+2); %1 x nv+2
d=floor(nx/2); % center point
if availGPU>0 && nx> 400 && gpu.AvailableMemory < 8e9
    device='gpuLowMem';
elseif nx<=100
    device='cpu';
end
%% Coefficients
disp('Compute Coefficients')
ticCoeffs=tic;
[A,B,Dx,Dv,Dxx,Dvv,H,Fx,Fv,Gxx,Gxv,Gvv,Sx,Sv,Phi]=...
    coefficients(xgrid,vgrid,h,fx,fv,gxx,gxv,gvv,sx,sv,phi);
ctimeCoeffs=toc(ticCoeffs);
fprintf('Elapsed time for coefficients is %3.3g seconds.\n',ctimeCoeffs)
% fprintf('Spectal norms: ||A||_2 = %3.3e, ||B||_2 = %3.3e,',norm(full(A),2),norm(full(B),2))
disp('Plot Sparsity Patterns')
if makePlots
    patterns=sparsityPattern(A,B,10000);
else
    patterns=[];
end
%% Switch to multithreading environment
% delete(gcp('nocreate'));parpool('threads');
%% 
phiEulerRef=repmat(Phi,1,1,1,M);
phiEuler=repmat(Phi,1,1,1,M);
phiMagnus1=repmat(Phi(:),1,M);
phiMagnus2=repmat(Phi(:),1,M);
phiMagnus3=repmat(Phi(:),1,M);
ctimeEulerRefTotal=0;
ctimeEulerTotal=0;
ctimeMagnus1Total=0;
ctimeMagnus2Total=0;
ctimeMagnus3Total=0;
ctimeExactTotal=0;
errDisExactEulerRef={};
errDisExactEuler={};
errDisExactMagnus1={};
errDisExactMagnus2={};
errDisExactMagnus3={};
errAbsExactEulerRef={};
errAbsExactEuler={};
errAbsExactMagnus1={};
errAbsExactMagnus2={};
errAbsExactMagnus3={};
errAbsAvgExactEulerRef={};
errAbsAvgExactEuler={};
errAbsAvgExactMagnus1={};
errAbsAvgExactMagnus2={};
errAbsAvgExactMagnus3={};

kk=1;
dWeulerRefCell=cell(floor(T/dT),1);
dWeulerCell=cell(floor(T/dT),1);
WmagnusCell=cell(floor(T/dT),1);
WexactCell=cell(floor(T/dT),1);
tiERMCell=cell(floor(T/dT),1);
tiERECell=cell(floor(T/dT),1);
tiEMCell=cell(floor(T/dT),1);
disp('Generate Brownian motion')
ticBM=tic;
for tk=dT:dT:T
[dWeulerRef,dWeuler,WeulerRef,Wmagnus,tiERM,tiERE,tiEM]=brownianMotion(dT,NeulerRef,Neuler,Nmagnus,M);
dWeulerRefCell{kk}=dWeulerRef;
dWeulerCell{kk}=dWeuler;
WmagnusCell{kk}=Wmagnus;
WexactCell{kk}=WeulerRef;
tiERMCell{kk}=tiERM;
tiERECell{kk}=tiERE;
tiEMCell{kk}=tiEM;
kk=kk+1;
end
ctimeBM=toc(ticBM);
fprintf('Elapsed time for BM is %3.3g seconds.\n',ctimeBM)
kk=1;
for tk=dT:dT:T
%% Brownian motion
disp('Select Brownian motion')
dWeulerRef=dWeulerRefCell{kk};
dWeuler=dWeulerCell{kk};
Wmagnus=WmagnusCell{kk};
tiERM=tiERMCell{kk};
tiERE=tiERECell{kk};
tiEM=tiEMCell{kk};

%% Time storage
%     tnMagnus=1:1:Nmagnus;
%     tnMagnus=[1,floor(Nmagnus/2),Nmagnus];
tnMagnus=Nmagnus;
tnEulerRef=tiERM(tnMagnus);
tnEuler=tiEM(tnMagnus);
%% Calculate reference Euler-Maruyama
if useEulerRef
    tEulerRef=linspace(0,dT,NeulerRef);
    disp('Compute reference Euler')
    ticEulerRef=tic;
    utEulerRef=eulerConst(tnEulerRef,dT,NeulerRef,dWeulerRef,...
                          Dx,Dv,Dxx,Dvv,H,Fx,Fv,Gxx,Gxv,Gvv,Sx,Sv,...
                          phiEulerRef,'device','cpu','mode','full');
    ctimeEulerRef=toc(ticEulerRef);
    ctimeEulerRefTotal=ctimeEulerRefTotal+ctimeEulerRef;
    fprintf('Elapsed time for Euler ref is %3.3g seconds.\n',ctimeEulerRef)
    phiEulerRef=utEulerRef(:,:,end,:);
    % Check if Euler scheme is diverging
    if any(isnan(utEulerRef),'all') || any(utEulerRef(floor(nx/2),floor(nv/2),:,:)>1e10,'all')
        disp('Diverging euler scheme')
    end
end
%% Calculate Euler-Maruyama
if useEuler
    tEuler=linspace(0,dT,Neuler);
    disp('Compute Euler')
    ticEuler=tic;
    utEuler=eulerConst(tnEuler,dT,Neuler,dWeuler,Dx,Dv,Dxx,Dvv,H,Fx,Fv,Gxx,Gxv,Gvv,Sx,Sv,phiEuler,'device','cpu','mode','full');
    ctimeEuler=toc(ticEuler);
    ctimeEulerTotal=ctimeEulerTotal+ctimeEuler;
    fprintf('Elapsed time for Euler is %3.3g seconds.\n',ctimeEuler)
    phiEuler=utEuler(:,:,end,:);
    % Check if Euler scheme is diverging
    if any(isnan(utEuler),'all') || any(utEuler(floor(nx/2),floor(nv/2),:,:)>1e10,'all')
    %     error('Diverging euler scheme')
        disp('Diverging euler scheme')
    end
end
%% Calculate Magnus
tMagnus=linspace(0,dT,Nmagnus);

if useMagnus1
    disp('Compute Magnus 1')
    ticMagnus1=tic;
    utMagnus1=magnusConst(tMagnus,Wmagnus,A,B,phiMagnus1,nx,nv,tnMagnus,1,'device',device);
    ctimeMagnus1=toc(ticMagnus1);
    fprintf('Elapsed time for Magnus order 1 is %3.3g seconds.\n',ctimeMagnus1)
    if any(isnan(utMagnus1),'all') || any(utMagnus1(floor(nx/2),floor(nv/2),:,:)>1e10,'all')
        disp('Diverging Magnus 1')
    end
    ctimeMagnus1Total=ctimeMagnus1Total+ctimeMagnus1;
    phiMagnus1=utMagnus1(:,:,end,:);
    phiMagnus1=reshape(phiMagnus1,[],M);
    if availGPU>0 && ~strcmp('device','cpu')
        clearGPU=parfevalOnAll(@gpuDevice,0,[]);
        wait(clearGPU)
    end
end
if useMagnus2
    disp('Compute Magnus 2')
    ticMagnus2=tic;
    utMagnus2=magnusConst(tMagnus,Wmagnus,A,B,phiMagnus2,nx,nv,tnMagnus,2,'device',device);
    ctimeMagnus2=toc(ticMagnus2);
    fprintf('Elapsed time for Magnus order 2 is %3.3g seconds.\n',ctimeMagnus2)
    if any(isnan(utMagnus2),'all') || any(utMagnus2(floor(nx/2),floor(nv/2),:,:)>1e10,'all')
        disp('Diverging Magnus 2')
    end
    ctimeMagnus2Total=ctimeMagnus2Total+ctimeMagnus2;
    phiMagnus2=utMagnus2(:,:,end,:);
    phiMagnus2=reshape(phiMagnus2,[],M);
    if availGPU>0 && ~strcmp('device','cpu')
        clearGPU=parfevalOnAll(@gpuDevice,0,[]);
        wait(clearGPU)
    end
end
%%
if useMagnus3
    disp('Compute Magnus 3')
    ticMagnus3=tic;
    utMagnus3=magnusConst(tMagnus,Wmagnus,A,B,phiMagnus3,nx,nv,tnMagnus,3,'device',device);
    ctimeMagnus3=toc(ticMagnus3);
    fprintf('Elapsed time for Magnus order 3 is %3.3g seconds.\n',ctimeMagnus3)
    if any(isnan(utMagnus3),'all') || any(utMagnus3(floor(nx/2),floor(nv/2),:,:)>1e10,'all')
        disp('Diverging Magnus 3')
    end
    ctimeMagnus3Total=ctimeMagnus3Total+ctimeMagnus3;
    phiMagnus3=utMagnus3(:,:,end,:);
    phiMagnus3=reshape(phiMagnus3,[],M);
    if availGPU>0 && ~strcmp('device','cpu')
        clearGPU=parfevalOnAll(@gpuDevice,0,[]);
        wait(clearGPU)
    end
end
%% Calculate Exact
if useExact && hasExact
    tExact=linspace(0,dT,NeulerRef);
    disp('Compute Exact')
    if kk==1
        tTemp=tExact;
        Wtemp=WexactCell{kk};
    else
        tTemp=cat(2,tTemp,tTemp(end)+tExact(2:end));
        Wtemp=cat(3,Wtemp,Wtemp(1,1,end,:)+WexactCell{kk}(1,1,2:end,:));
    end
    Ti=length(tTemp);
    ticExact=tic;
    utExact=exact(tk,Ti,tTemp,xgrid(2:end-1),vgrid(2:end-1),a,sigma,Wtemp);
    ctimeExact=toc(ticExact);
    ctimeExactTotal=ctimeExactTotal+ctimeExact;
    fprintf('Elapsed time for Exact is %3.3g seconds.\n',ctimeExact)
end
%% Calculate Errors
disp('Compute errors')
if useExact && hasExact
    disp('Reference Method exact')
    if useEulerRef
        [errAbsExactEulerRef{kk},errAbsAvgExactEulerRef{kk}]=...
            absError(utExact,utEulerRef,kappas);
        fprintf('\t Mean abs error Exact and Euler Ref at time %1.3g at center: %1.3e\n',...
                tk,errAbsExactEulerRef{kk}(d,d));
        [errDisExactEulerRef{kk},xExactEulerRef{kk},fExactEulerRef{kk},...
            nBlowupsExactEulerRef{kk}]= ...
            errorDis(utExact,utEulerRef,kappas,...
                    'blowup',blowup,...
                    'normType',errorDisType);
        for iKappa=1:1:length(kappas)
            if nBlowupsExactEulerRef{kk}{iKappa}>-1
                fprintf('\t\t Number of blowups for Euler Ref with kappa=%d: %d\n',...
                        kappas(iKappa),nBlowupsExactEulerRef{kk}{iKappa})
            end
        end
    end
    if useEuler
        [errAbsExactEuler{kk},errAbsAvgExactEuler{kk}]=...
            absError(utExact,utEuler,kappas);
        fprintf('\t Mean abs error Exact and Euler at time %1.3g at center: %1.3e\n',...
                    tk,errAbsExactEuler{kk}(d,d));
        [errDisExactEuler{kk},xExactEuler{kk},fExactEuler{kk},...
            nBlowupsExactEuler{kk}]= ...
            errorDis(utExact,utEuler,kappas,...
                    'blowup',blowup,...
                    'normType',errorDisType);
        for iKappa=1:1:length(kappas)
            if nBlowupsExactEuler{kk}{iKappa}>-1
                fprintf('\t\t Number of blowups for Euler with kappa=%d: %d\n',...
                        kappas(iKappa),nBlowupsExactEuler{kk}{iKappa})
            end
        end
    end
    if useMagnus1
        [errAbsExactMagnus1{kk},errAbsAvgExactMagnus1{kk}]=...
            absError(utExact,utMagnus1,kappas);
        fprintf('\t Mean abs error Exact and Magnus 1 at time %1.3g at center: %1.3e\n',...
                    tk,errAbsExactMagnus1{kk}(d,d));
        [errDisExactMagnus1{kk},xExactMagnus1{kk},fExactMagnus1{kk},...
            nBlowupsExactMagnus1{kk}]= ...
            errorDis(utExact,utMagnus1,kappas,...
                    'blowup',blowup,...
                    'normType',errorDisType);
        for iKappa=1:1:length(kappas)
            if nBlowupsExactMagnus1{kk}{iKappa}>-1
                fprintf('\t\t Number of blowups for Magnus 1 with kappa=%d: %d\n',...
                        kappas(iKappa),nBlowupsExactMagnus1{kk}{iKappa})
            end
        end
    end
    if useMagnus2
        [errAbsExactMagnus2{kk},errAbsAvgExactMagnus2{kk}]=...
            absError(utExact,utMagnus2,kappas);
        fprintf('\t Mean abs error Exact and Magnus 2 at time %1.3g at center: %1.3e\n',...
                    tk,errAbsExactMagnus2{kk}(d,d));
        [errDisExactMagnus2{kk},xExactMagnus2{kk},fExactMagnus2{kk},...
            nBlowupsExactMagnus2{kk}]= ...
            errorDis(utExact,utMagnus2,kappas,...
                    'blowup',blowup,...
                    'normType',errorDisType);
        for iKappa=1:1:length(kappas)
            if nBlowupsExactMagnus2{kk}{iKappa}>-1
                fprintf('\t\t Number of blowups for Magnus 2 with kappa=%d: %d\n',...
                        kappas(iKappa),nBlowupsExactMagnus2{kk}{iKappa})
            end
        end
    end
    if useMagnus3
        [errAbsExactMagnus3{kk},errAbsAvgExactMagnus3{kk}]=...
            absError(utExact,utMagnus3,kappas);
        fprintf('\t Mean abs error Exact and Magnus 3 at time %1.3g at center: %1.3e\n',...
                    tk,errAbsExactMagnus3{kk}(d,d));
        [errDisExactMagnus3{kk},xExactMagnus3{kk},fExactMagnus3{kk},...
            nBlowupsExactMagnus3{kk}]= ...
            errorDis(utExact,utMagnus3,kappas,...
                    'blowup',blowup,...
                    'normType',errorDisType);
        for iKappa=1:1:length(kappas)
            if nBlowupsExactMagnus3{kk}{iKappa}>-1
                fprintf('\t\t Number of blowups for Magnus 3 with kappa=%d: %d\n',...
                        kappas(iKappa),nBlowupsExactMagnus3{kk}{iKappa})
            end
        end
    end
end
if useEulerRef
    disp('Reference Method EulerRef')
    if useEuler
        [errAbsEulerRefEuler{kk},errAbsAvgEulerRefEuler{kk}]=...
            absError(utEulerRef,utEuler,kappas);
        fprintf('\t Mean abs error Euler Ref and Euler at time %1.3g at center: %1.3e\n',...
                        tk,errAbsEulerRefEuler{kk}(d,d));
        [errDisEulerRefEuler{kk},xEulerRefEuler{kk},fEulerRefEuler{kk},...
            nBlowupsEulerRefEuler{kk}]= ...
            errorDis(utEulerRef,utEuler,kappas,...
                    'blowup',blowup,...
                    'normType',errorDisType);
        for iKappa=1:1:length(kappas)
            if nBlowupsEulerRefEuler{kk}{iKappa}>-1
                fprintf('\t\t Number of blowups for Euler with kappa=%d: %d\n',...
                        kappas(iKappa),nBlowupsEulerRefEuler{kk}{iKappa})
            end
        end
    end
    if useMagnus1
        [errAbsEulerRefMagnus1{kk},errAbsAvgEulerRefMagnus1{kk}]=...
            absError(utEulerRef,utMagnus1,kappas);
        fprintf('\t Mean abs error Euler Ref and Magnus 1 at time %1.3g at center: %1.3e\n',...
                        tk,errAbsEulerRefMagnus1{kk}(d,d));
        [errDisEulerRefMagnus1{kk},xEulerRefMagnus1{kk},fEulerRefMagnus1{kk},...
            nBlowupsEulerRefMagnus1{kk}]= ...
            errorDis(utEulerRef,utMagnus1,kappas,...
                    'blowup',blowup,...
                    'normType',errorDisType);
        for iKappa=1:1:length(kappas)
            if nBlowupsEulerRefMagnus1{kk}{iKappa}>-1
                fprintf('\t\t Number of blowups for Magnus 1 with kappa=%d: %d\n',...
                        kappas(iKappa),nBlowupsEulerRefMagnus1{kk}{iKappa})
            end
        end
    end
    if useMagnus2
        [errAbsEulerRefMagnus2{kk},errAbsAvgEulerRefMagnus2{kk}]=...
            absError(utEulerRef,utMagnus2,kappas);
        fprintf('\t Mean abs error Euler Ref and Magnus 2 at time %1.3g at center: %1.3e\n',...
                        tk,errAbsEulerRefMagnus2{kk}(d,d));
        [errDisEulerRefMagnus2{kk},xEulerRefMagnus2{kk},fEulerRefMagnus2{kk},...
            nBlowupsEulerRefMagnus2{kk}]= ...
            errorDis(utEulerRef,utMagnus2,kappas,...
                    'blowup',blowup,...
                    'normType',errorDisType);
        for iKappa=1:1:length(kappas)
            if nBlowupsEulerRefMagnus2{kk}{iKappa}>-1
                fprintf('\t\t Number of blowups for Magnus 2 with kappa=%d: %d\n',...
                        kappas(iKappa),nBlowupsEulerRefMagnus2{kk}{iKappa})
            end
        end
    end
    if useMagnus3
        [errAbsEulerRefMagnus3{kk},errAbsAvgEulerRefMagnus3{kk}]=...
            absError(utEulerRef,utMagnus3,kappas);
        fprintf('\t Mean abs error Euler Ref and Magnus 3 at time %1.3g at center: %1.3e\n',...
                        tk,errAbsEulerRefMagnus3{kk}(d,d));
        [errDisEulerRefMagnus3{kk},xEulerRefMagnus3{kk},fEulerRefMagnus3{kk},...
            nBlowupsEulerRefMagnus3{kk}]= ...
            errorDis(utEulerRef,utMagnus3,kappas,...
                    'blowup',blowup,...
                    'normType',errorDisType);
        for iKappa=1:1:length(kappas)
            if nBlowupsEulerRefMagnus3{kk}{iKappa}>-1
                fprintf('\t\t Number of blowups for Magnus 3 with kappa=%d: %d\n',...
                        kappas(iKappa),nBlowupsEulerRefMagnus3{kk}{iKappa})
            end
        end
    end
end
kk=kk+1;
end
%% FileName
fileName=sprintf('%s_T%1.3f_d%d_NER%d_NE_%d_NM_%d_M%d_%d%d%d%d%d%d',...
    methodName,T,nx,NeulerRef,Neuler,floor(T/dT)+1,M,...
    hasExact,useEulerRef,useEuler,useMagnus1,useMagnus2,useMagnus3);
%% Total computational times
fprintf('Total Elapsed time for Euler ref is %3.3g seconds.\n',ctimeEulerRefTotal)
fprintf('Total Elapsed time for Euler is %3.3g seconds.\n',ctimeEulerTotal)
fprintf('Total Elapsed time for Magnus order 1 is %3.3g seconds.\n',ctimeMagnus1Total)
fprintf('Total Elapsed time for Magnus order 2 is %3.3g seconds.\n',ctimeMagnus2Total)
fprintf('Total Elapsed time for Magnus order 3 is %3.3g seconds.\n',ctimeMagnus3Total)
fprintf('Total Elapsed time for Exact is %3.3g seconds.\n',ctimeExactTotal)
%% Error plots 
tEval=dT:dT:T;
showTitle=true;
showLines=true;
absErrorExactFig=[];
absErrorExactVid={};
errorDisExactFig=[];
errorDisExactVid={};
if useExact && hasExact && makePlots
    if useEulerRef
        tempFig=plotAbsErrors(errAbsExactEulerRef,errAbsAvgExactEulerRef,...
                          xgrid(2:end-1),vgrid(2:end-1),kappas,tEval,...
                          'showTitle',showTitle,...
                          'showLines',showLines, ...
                          'backgroundColor',backgroundColor,...
                          'textColor',textColor);
        tempFig2=plotErrorDis(errDisExactEulerRef,xExactEulerRef,...
                              fExactEulerRef,nBlowupsExactEulerRef,...
                              kappas,tEval,...
                              'showTitle',showTitle,...
                              'percentile',percentile, ...
                              'backgroundColor',backgroundColor,...
                              'textColor',textColor);
        if makeErrorGif
            v=videoAbsErrors([dirAbsErrorGif,'/',fileName],'ExactEulerRef',tempFig, ...
                          'backgroundColor',backgroundColor,...
                          'textColor',textColor);
            absErrorExactVid{end+1}=v;
            v=videoErrorDis([dirErrorDisGif,'/',fileName],'ExactEulerRef',tempFig2, ...
                          'backgroundColor',backgroundColor,...
                          'textColor',textColor);
            errorDisExactVid=cat(2,errorDisExactVid,v);
        end
        absErrorExactFig=cat(2,absErrorExactFig,...
            [tempFig{1},tempFig{floor(length(tempFig)/2)},tempFig{end}]);
        tempFig21=[tempFig2{1,:},tempFig2{floor(length(tempFig)/2),:},tempFig2{end,:}];
        errorDisExactFig=cat(2,errorDisExactFig,tempFig21);
        close(setdiff(cat(2,tempFig{:}),absErrorExactFig))
        close(setdiff(cat(2,tempFig2{:}),errorDisExactFig))
    end
    if useMagnus1
        tempFig=plotAbsErrors(errAbsExactMagnus1,errAbsAvgExactMagnus1,...
                          xgrid(2:end-1),vgrid(2:end-1),kappas,tEval,...
                          'showTitle',showTitle,...
                          'showLines',showLines, ...
                          'backgroundColor',backgroundColor,...
                          'textColor',textColor);
        tempFig2=plotErrorDis(errDisExactMagnus1,xExactMagnus1,...
                              fExactMagnus1,nBlowupsExactMagnus1,...
                              kappas,tEval,...
                              'showTitle',showTitle,...
                              'percentile',percentile, ...
                          'backgroundColor',backgroundColor,...
                          'textColor',textColor);
        if makeErrorGif
            v=videoAbsErrors([dirAbsErrorGif,'/',fileName],'ExactMagnus1',tempFig, ...
                          'backgroundColor',backgroundColor,...
                          'textColor',textColor);
            absErrorExactVid{end+1}=v;
            v=videoErrorDis([dirErrorDisGif,'/',fileName],'ExactMagnus1',tempFig2, ...
                          'backgroundColor',backgroundColor,...
                          'textColor',textColor);
            errorDisExactVid=cat(2,errorDisExactVid,v);
        end
        absErrorExactFig=cat(2,absErrorExactFig,...
            [tempFig{1},tempFig{floor(length(tempFig)/2)},tempFig{end}]);
        tempFig21=[tempFig2{1,:},tempFig2{floor(length(tempFig)/2),:},tempFig2{end,:}];
        errorDisExactFig=cat(2,errorDisExactFig,tempFig21);
        close(setdiff(cat(2,tempFig{:}),absErrorExactFig))
        close(setdiff(cat(2,tempFig2{:}),errorDisExactFig))
    end
    if useMagnus2
        tempFig=plotAbsErrors(errAbsExactMagnus2,errAbsAvgExactMagnus2,...
                          xgrid(2:end-1),vgrid(2:end-1),kappas,tEval,...
                          'showTitle',showTitle,...
                          'showLines',showLines, ...
                          'backgroundColor',backgroundColor,...
                          'textColor',textColor);
        tempFig2=plotErrorDis(errDisExactMagnus2,xExactMagnus2,...
                              fExactMagnus2,nBlowupsExactMagnus2,...
                              kappas,tEval,...
                              'showTitle',showTitle,...
                              'percentile',percentile, ...
                          'backgroundColor',backgroundColor,...
                          'textColor',textColor);
        if makeErrorGif
            v=videoAbsErrors([dirAbsErrorGif,'/',fileName],'ExactMagnus2',tempFig, ...
                          'backgroundColor',backgroundColor,...
                          'textColor',textColor);
            absErrorExactVid{end+1}=v;
            v=videoErrorDis([dirErrorDisGif,'/',fileName],'ExactMagnus2',tempFig2, ...
                          'backgroundColor',backgroundColor,...
                          'textColor',textColor);
            errorDisExactVid=cat(2,errorDisExactVid,v);
        end
        absErrorExactFig=cat(2,absErrorExactFig,...
            [tempFig{1},tempFig{floor(length(tempFig)/2)},tempFig{end}]);
        tempFig21=[tempFig2{1,:},tempFig2{floor(length(tempFig)/2),:},tempFig2{end,:}];
        errorDisExactFig=cat(2,errorDisExactFig,tempFig21);
        close(setdiff(cat(2,tempFig{:}),absErrorExactFig))
        close(setdiff(cat(2,tempFig2{:}),errorDisExactFig))
    end
    if useMagnus3
        tempFig=plotAbsErrors(errAbsExactMagnus3,errAbsAvgExactMagnus3,...
                          xgrid(2:end-1),vgrid(2:end-1),kappas,tEval,...
                          'showTitle',showTitle,...
                          'showLines',showLines, ...
                          'backgroundColor',backgroundColor,...
                          'textColor',textColor);
        tempFig2=plotErrorDis(errDisExactMagnus3,xExactMagnus3,...
                              fExactMagnus3,nBlowupsExactMagnus3,...
                              kappas,tEval,...
                              'showTitle',showTitle,...
                              'percentile',percentile, ...
                          'backgroundColor',backgroundColor,...
                          'textColor',textColor);
        if makeErrorGif
            v=videoAbsErrors([dirAbsErrorGif,'/',fileName],'ExactMagnus3',tempFig, ...
                          'backgroundColor',backgroundColor,...
                          'textColor',textColor);
            absErrorExactVid{end+1}=v;
            v=videoErrorDis([dirErrorDisGif,'/',fileName],'ExactMagnus3',tempFig2, ...
                          'backgroundColor',backgroundColor,...
                          'textColor',textColor);
            errorDisExactVid=cat(2,errorDisExactVid,v);
        end
        absErrorExactFig=cat(2,absErrorExactFig,...
            [tempFig{1},tempFig{floor(length(tempFig)/2)},tempFig{end}]);
        tempFig21=[tempFig2{1,:},tempFig2{floor(length(tempFig)/2),:},tempFig2{end,:}];
        errorDisExactFig=cat(2,errorDisExactFig,tempFig21);
        close(setdiff(cat(2,tempFig{:}),absErrorExactFig))
        close(setdiff(cat(2,tempFig2{:}),errorDisExactFig))
    end
end
%%
absErrorEulerRefFig=[];
absErrorEulerRefVid={};
errorDisEulerRefFig=[];
errorDisEulerRefVid={};
if useEulerRef && ~(useExact && hasExact) && makePlots
    if useEuler
        tempFig=plotAbsErrors(errAbsEulerRefEuler,errAbsAvgEulerRefEuler,...
                          xgrid(2:end-1),vgrid(2:end-1),kappas,tEval,...
                          'showTitle',showTitle,...
                          'showLines',showLines, ...
                          'backgroundColor',backgroundColor,...
                          'textColor',textColor);
        tempFig2=plotErrorDis(errDisEulerRefEuler,xEulerRefEuler,...
                              fEulerRefEuler,nBlowupsEulerRefEuler,...
                              kappas,tEval,...
                              'showTitle',showTitle,...
                              'percentile',percentile, ...
                          'backgroundColor',backgroundColor,...
                          'textColor',textColor);
        if makeErrorGif
            v=videoAbsErrors([dirAbsErrorGif,'/',fileName],'EulerRefEuler',tempFig, ...
                          'backgroundColor',backgroundColor,...
                          'textColor',textColor);
            absErrorEulerRefVid{end+1}=v;
            v=videoErrorDis([dirErrorDisGif,'/',fileName],'EulerRefEuler',tempFig2, ...
                          'backgroundColor',backgroundColor,...
                          'textColor',textColor);
            errorDisEulerRefVid=cat(2,errorDisEulerRefVid,v);
        end
        absErrorEulerRefFig=cat(2,absErrorEulerRefFig,...
            [tempFig{1},tempFig{floor(length(tempFig)/2)},tempFig{end}]);
        tempFig21=[tempFig2{1,:},tempFig2{floor(length(tempFig)/2),:},tempFig2{end,:}];
        errorDisEulerRefFig=cat(2,errorDisEulerRefFig,tempFig21);
        close(setdiff(cat(2,tempFig{:}),absErrorEulerRefFig))
        close(setdiff(cat(2,tempFig2{:}),errorDisEulerRefFig))
    end
    if useMagnus1
        tempFig=plotAbsErrors(errAbsEulerRefMagnus1,errAbsAvgEulerRefMagnus1,...
                          xgrid(2:end-1),vgrid(2:end-1),kappas,tEval,...
                          'showTitle',showTitle,...
                          'showLines',showLines, ...
                          'backgroundColor',backgroundColor,...
                          'textColor',textColor);
        tempFig2=plotErrorDis(errDisEulerRefMagnus1,xEulerRefMagnus1,...
                              fEulerRefMagnus1,nBlowupsEulerRefMagnus1,...
                              kappas,tEval,...
                              'showTitle',showTitle,...
                              'percentile',percentile, ...
                          'backgroundColor',backgroundColor,...
                          'textColor',textColor);
        if makeErrorGif
            v=videoAbsErrors([dirAbsErrorGif,'/',fileName],'EulerRefMagnus1',tempFig, ...
                          'backgroundColor',backgroundColor,...
                          'textColor',textColor);
            absErrorEulerRefVid{end+1}=v;
            v=videoErrorDis([dirErrorDisGif,'/',fileName],'EulerRefMagnus1',tempFig2, ...
                          'backgroundColor',backgroundColor,...
                          'textColor',textColor);
            errorDisEulerRefVid=cat(2,errorDisEulerRefVid,v);
        end
        absErrorEulerRefFig=cat(2,absErrorEulerRefFig,...
            [tempFig{1},tempFig{floor(length(tempFig)/2)},tempFig{end}]);
        tempFig21=[tempFig2{1,:},tempFig2{floor(length(tempFig)/2),:},tempFig2{end,:}];
        errorDisEulerRefFig=cat(2,errorDisEulerRefFig,tempFig21);
        close(setdiff(cat(2,tempFig{:}),absErrorEulerRefFig))
        close(setdiff(cat(2,tempFig2{:}),errorDisEulerRefFig))
    end
    if useMagnus2
        tempFig=plotAbsErrors(errAbsEulerRefMagnus2,errAbsAvgEulerRefMagnus2,...
                          xgrid(2:end-1),vgrid(2:end-1),kappas,tEval,...
                          'showTitle',showTitle,...
                          'showLines',showLines, ...
                          'backgroundColor',backgroundColor,...
                          'textColor',textColor);
        tempFig2=plotErrorDis(errDisEulerRefMagnus2,xEulerRefMagnus2,...
                              fEulerRefMagnus2,nBlowupsEulerRefMagnus2,...
                              kappas,tEval,...
                              'showTitle',showTitle,...
                              'percentile',percentile, ...
                          'backgroundColor',backgroundColor,...
                          'textColor',textColor);
        if makeErrorGif
            v=videoAbsErrors([dirAbsErrorGif,'/',fileName],'EulerRefMagnus2',tempFig, ...
                          'backgroundColor',backgroundColor,...
                          'textColor',textColor);
            absErrorEulerRefVid{end+1}=v;
            v=videoErrorDis([dirErrorDisGif,'/',fileName],'EulerRefMagnus2',tempFig2, ...
                          'backgroundColor',backgroundColor,...
                          'textColor',textColor);
            errorDisEulerRefVid=cat(2,errorDisEulerRefVid,v);
        end
        absErrorEulerRefFig=cat(2,absErrorEulerRefFig,...
            [tempFig{1},tempFig{floor(length(tempFig)/2)},tempFig{end}]);
        tempFig21=[tempFig2{1,:},tempFig2{floor(length(tempFig)/2),:},tempFig2{end,:}];
        errorDisEulerRefFig=cat(2,errorDisEulerRefFig,tempFig21);
        close(setdiff(cat(2,tempFig{:}),absErrorEulerRefFig))
        close(setdiff(cat(2,tempFig2{:}),errorDisEulerRefFig))
    end
    if useMagnus3
        tempFig=plotAbsErrors(errAbsEulerRefMagnus3,errAbsAvgEulerRefMagnus3,...
                          xgrid(2:end-1),vgrid(2:end-1),kappas,tEval,...
                          'showTitle',showTitle,...
                          'showLines',showLines, ...
                          'backgroundColor',backgroundColor,...
                          'textColor',textColor);
        tempFig2=plotErrorDis(errDisEulerRefMagnus3,xEulerRefMagnus3,...
                              fEulerRefMagnus3,nBlowupsEulerRefMagnus3,...
                              kappas,tEval,...
                              'showTitle',showTitle,...
                              'percentile',percentile, ...
                          'backgroundColor',backgroundColor,...
                          'textColor',textColor);
        if makeErrorGif
            v=videoAbsErrors([dirAbsErrorGif,'/',fileName],'EulerRefMagnus3',tempFig, ...
                          'backgroundColor',backgroundColor,...
                          'textColor',textColor);
            absErrorEulerRefVid{end+1}=v;
            v=videoErrorDis([dirErrorDisGif,'/',fileName],'EulerRefMagnus3',tempFig2, ...
                          'backgroundColor',backgroundColor,...
                          'textColor',textColor);
            errorDisEulerRefVid=cat(2,errorDisEulerRefVid,v);
        end
        absErrorEulerRefFig=cat(2,absErrorEulerRefFig,...
            [tempFig{1},tempFig{floor(length(tempFig)/2)},tempFig{end}]);
        tempFig21=[tempFig2{1,:},tempFig2{floor(length(tempFig)/2),:},tempFig2{end,:}];
        errorDisEulerRefFig=cat(2,errorDisEulerRefFig,tempFig21);
        close(setdiff(cat(2,tempFig{:}),absErrorEulerRefFig))
        close(setdiff(cat(2,tempFig2{:}),errorDisEulerRefFig))
    end
end
%% Output
disp('Create output')
output(fileName,methodName,xl,xu,nx,vl,vu,nv,T,tnEuler,tnMagnus,dtEulerRef,dtEuler,dT,M,...
                h,fx,fv,gxx,gxv,gvv,sx,sv,phi,...
                ctimeExactTotal,ctimeEulerRefTotal,ctimeEulerTotal,ctimeMagnus1Total,ctimeMagnus2Total,ctimeMagnus3Total,...
                xgrid(2:end-1),vgrid(2:end-1),kappas,...
                errDisExactEulerRef,errDisExactEuler,errDisExactMagnus1,errDisExactMagnus2,errDisExactMagnus3,...
                errDisEulerRefEuler,errDisEulerRefMagnus1,errDisEulerRefMagnus2,errDisEulerRefMagnus3,...
                [],[absErrorExactFig,errorDisExactFig],[absErrorEulerRefFig,errorDisEulerRefFig],patterns,...
                hasExact,useExact,useEulerRef,useEuler,useMagnus1,useMagnus2,useMagnus3,...
                tEval,absErrorExactVid,absErrorEulerRefVid,errorDisExactVid,errorDisEulerRefVid,...
                errAbsAvgExactEulerRef,errAbsAvgExactEuler,errAbsAvgExactMagnus1,errAbsAvgExactMagnus2,errAbsAvgExactMagnus3,...
                errAbsAvgEulerRefEuler,errAbsAvgEulerRefMagnus1,errAbsAvgEulerRefMagnus2,errAbsAvgEulerRefMagnus3,...
                'backgroundColor',backgroundColor);
outputMatlab(fileName,...
              T,dtEulerRef,dtEuler,dT,M,kappas,nx,nv,...
              ctimeEulerRefTotal,ctimeEulerTotal,ctimeMagnus1Total,ctimeMagnus2Total,ctimeMagnus3Total,...
              errDisExactEulerRef,errDisExactEuler,errDisExactMagnus1,errDisExactMagnus2,errDisExactMagnus3,...
              errDisEulerRefEuler,errDisEulerRefMagnus1,errDisEulerRefMagnus2,errDisEulerRefMagnus3,...
              errAbsAvgExactEulerRef,errAbsAvgExactEuler,errAbsAvgExactMagnus1,errAbsAvgExactMagnus2,errAbsAvgExactMagnus3,...
              errAbsAvgEulerRefEuler,errAbsAvgEulerRefMagnus1,errAbsAvgEulerRefMagnus2,errAbsAvgEulerRefMagnus3,...
              hasExact,useExact,useEulerRef,useEuler,useMagnus1,useMagnus2,useMagnus3)
disp('done')
end