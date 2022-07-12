clear all; close all; fclose('all'); rng(0);
availGPU=gpuDeviceCount("available");
if availGPU > 0
    device='gpu';
    gpu=gpuDevice(1);
else
    device='cpu';
end
delete(gcp('nocreate'));parpool('threads'); % multithreading
%% Select methods
useEulerRef = true;
useMagnus1 = true;
useMagnus2 = true;
useMagnus3 = true;
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
% [h,fx,fv,gxx,gxv,gvv,sx,sv,phi,hasExact,methodName]=...
%     selectCoefficients('exact1',[a,sigma]);
[h,fx,fv,gxx,gxv,gvv,sx,sv,phi,hasExact,methodName]=...
    selectCoefficients('Langevin1',[a,sigma]);
%% Video & Plot parameters
backgroundColor = 'w';
% backgroundColor = [53,54,58]./255;
textColor = 'k';
% textColor = [237,237,237]./255;

%% Simulations
M = 100;
%% Space grids
% Position
xl=-4;
xu=4;
nx=200;
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

%%
T=1;
dtMagnusLog=1e-4;
dtEulerRef=1e-4;
% for dT = 0.01:0.01:0.02
dTtemp = unique([T./(2:1:100),T./(2.^(1:10)),T./(10.^(1:3))]);
dTtemp = dTtemp(3:3:end-1);
for dT = dTtemp
    %% Initialize variables
    ctimeMagnus1Total=0;
    ctimeMagnus2Total=0;
    ctimeMagnus3Total=0;
    ctimeEulerRefTotal=0;
    errDisEulerRefMagnus1={};
    errDisEulerRefMagnus2={};
    errDisEulerRefMagnus3={};
    errAbsEulerRefMagnus1={};
    errAbsEulerRefMagnus2={};
    errAbsEulerRefMagnus3={};
    errAbsAvgEulerRefMagnus1={};
    errAbsAvgEulerRefMagnus2={};
    errAbsAvgEulerRefMagnus3={};
    xEulerRefMagnus1={};
    fEulerRefMagnus1={};
    nBlowupsEulerRefMagnus1={};
    xEulerRefMagnus2={};
    fEulerRefMagnus2={};
    nBlowupsEulerRefMagnus2={};
    xEulerRefMagnus3={};
    fEulerRefMagnus3={};
    nBlowupsEulerRefMagnus3={};
    
    NeulerRef = floor(dT/dtEulerRef)+1;
    Nmagnus = floor(dT/dtMagnusLog)+1;
    phiEulerRef=repmat(Phi,1,1,1,M);
    phiMagnus1=repmat(Phi(:),1,M);
    phiMagnus2=repmat(Phi(:),1,M);
    phiMagnus3=repmat(Phi(:),1,M);
    
    kk=1;
    WmagnusCell=cell(floor(T/dT),1);
    dWeulerRefCell=cell(floor(T/dT),1);
    tiERMCell=cell(floor(T/dT),1);
%     disp('Generate Brownian motion')
    for tk=dT:dT:T
%         [~,~,Wmagnus,~,~,~,~]=brownianMotion(dT,Nmagnus,1,1,M);
        [dWeulerRef,~,~,Wmagnus,tiERM,~,~]=brownianMotion(dT,NeulerRef,1,Nmagnus,M);
        WmagnusCell{kk}=Wmagnus;
        dWeulerRefCell{kk}=dWeulerRef;
        tiERMCell{kk}=tiERM;
        kk=kk+1;
    end
    kk=1;
    fprintf('Curr dT=%g\n',dT)
    for tk=dT:dT:T
        %% Brownian motion
        tnMagnus=Nmagnus;
        tnEulerRef=tiERM(tnMagnus);
        Wmagnus=WmagnusCell{kk};
        dWeulerRef=dWeulerRefCell{kk};
        %% Calculate reference Euler-Maruyama
        if useEulerRef
            tEulerRef=linspace(0,dT,NeulerRef);
%             disp('Compute reference Euler')
            ticEulerRef=tic;
            utEulerRef=eulerConst(tnEulerRef,dT,NeulerRef,dWeulerRef,...
                                  Dx,Dv,Dxx,Dvv,H,Fx,Fv,Gxx,Gxv,Gvv,Sx,Sv,...
                                  phiEulerRef,'device','cpu','mode','full');
            ctimeEulerRef=toc(ticEulerRef);
            ctimeEulerRefTotal=ctimeEulerRefTotal+ctimeEulerRef;
%             fprintf('Elapsed time for Euler ref is %3.3g seconds.\n',ctimeEulerRef)
            phiEulerRef=utEulerRef(:,:,end,:);
            % Check if Euler scheme is diverging
%             if any(isnan(utEulerRef),'all') || any(utEulerRef(floor(nx/2),floor(nv/2),:,:)>1e10,'all')
%                 disp('Diverging euler scheme')
%             end
        end
        %% Calculate Magnus
        tMagnus=linspace(0,dT,Nmagnus);
        
        if useMagnus1
%             disp('Compute Magnus 1')
            ticMagnus1=tic;
            utMagnus1=magnusConst(tMagnus,Wmagnus,A,B,phiMagnus1,nx,nv,tnMagnus,1,'device',device);
            ctimeMagnus1=toc(ticMagnus1);
%             fprintf('Elapsed time for Magnus order 1 is %3.3g seconds.\n',ctimeMagnus1)
%             if any(isnan(utMagnus1),'all') || any(utMagnus1(floor(nx/2),floor(nv/2),:,:)>1e10,'all')
%                 disp('Diverging Magnus 1')
%             end
            ctimeMagnus1Total=ctimeMagnus1Total+ctimeMagnus1;
            phiMagnus1=utMagnus1(:,:,end,:);
            phiMagnus1=reshape(phiMagnus1,[],M);
            if availGPU>0 && ~strcmp('device','cpu')
                clearGPU=parfevalOnAll(@gpuDevice,0,[]);
                wait(clearGPU)
            end
        end
        if useMagnus2
%             disp('Compute Magnus 2')
            ticMagnus2=tic;
            utMagnus2=magnusConst(tMagnus,Wmagnus,A,B,phiMagnus2,nx,nv,tnMagnus,2,'device',device);
            ctimeMagnus2=toc(ticMagnus2);
%             fprintf('Elapsed time for Magnus order 2 is %3.3g seconds.\n',ctimeMagnus2)
%             if any(isnan(utMagnus2),'all') || any(utMagnus2(floor(nx/2),floor(nv/2),:,:)>1e10,'all')
%                 disp('Diverging Magnus 2')
%             end
            ctimeMagnus2Total=ctimeMagnus2Total+ctimeMagnus2;
            phiMagnus2=utMagnus2(:,:,end,:);
            phiMagnus2=reshape(phiMagnus2,[],M);
            if availGPU>0 && ~strcmp('device','cpu')
                clearGPU=parfevalOnAll(@gpuDevice,0,[]);
                wait(clearGPU)
            end
        end
        if useMagnus3
%             disp('Compute Magnus 3')
            ticMagnus3=tic;
            utMagnus3=magnusConst(tMagnus,Wmagnus,A,B,phiMagnus3,nx,nv,tnMagnus,3,'device',device);
            ctimeMagnus3=toc(ticMagnus3);
%             fprintf('Elapsed time for Magnus order 3 is %3.3g seconds.\n',ctimeMagnus3)
%             if any(isnan(utMagnus3),'all') || any(utMagnus3(floor(nx/2),floor(nv/2),:,:)>1e10,'all')
%                 disp('Diverging Magnus 3')
%             end
            ctimeMagnus3Total=ctimeMagnus3Total+ctimeMagnus3;
            phiMagnus3=utMagnus3(:,:,end,:);
            phiMagnus3=reshape(phiMagnus3,[],M);
            if availGPU>0 && ~strcmp('device','cpu')
                clearGPU=parfevalOnAll(@gpuDevice,0,[]);
                wait(clearGPU)
            end
        end
        %% Errors
        if useMagnus1
            [errAbsEulerRefMagnus1{kk},errAbsAvgEulerRefMagnus1{kk}]=...
                absError(utEulerRef,utMagnus1,kappas);
%             fprintf('\t Mean abs error Euler Ref and Magnus 1 at time %1.3g at center: %1.3e\n',...
%                             tk,errAbsEulerRefMagnus1{kk}(d,d));
            [errDisEulerRefMagnus1{kk},xEulerRefMagnus1{kk},fEulerRefMagnus1{kk},...
                nBlowupsEulerRefMagnus1{kk}]= ...
                errorDis(utEulerRef,utMagnus1,kappas,...
                        'blowup',blowup,...
                        'normType',errorDisType);
%             for iKappa=1:1:length(kappas)
%                 if nBlowupsEulerRefMagnus1{kk}{iKappa}>-1
%                     fprintf('\t\t Number of blowups for Magnus 1 with kappa=%d: %d\n',...
%                             kappas(iKappa),nBlowupsEulerRefMagnus1{kk}{iKappa})
%                 end
%             end
        end
        if useMagnus2
            [errAbsEulerRefMagnus2{kk},errAbsAvgEulerRefMagnus2{kk}]=...
                absError(utEulerRef,utMagnus2,kappas);
%             fprintf('\t Mean abs error Euler Ref and Magnus 2 at time %1.3g at center: %1.3e\n',...
%                             tk,errAbsEulerRefMagnus2{kk}(d,d));
            [errDisEulerRefMagnus2{kk},xEulerRefMagnus2{kk},fEulerRefMagnus2{kk},...
                nBlowupsEulerRefMagnus2{kk}]= ...
                errorDis(utEulerRef,utMagnus2,kappas,...
                        'blowup',blowup,...
                        'normType',errorDisType);
%             for iKappa=1:1:length(kappas)
%                 if nBlowupsEulerRefMagnus2{kk}{iKappa}>-1
%                     fprintf('\t\t Number of blowups for Magnus 2 with kappa=%d: %d\n',...
%                             kappas(iKappa),nBlowupsEulerRefMagnus2{kk}{iKappa})
%                 end
%             end
        end
        if useMagnus3
            [errAbsEulerRefMagnus3{kk},errAbsAvgEulerRefMagnus3{kk}]=...
                absError(utEulerRef,utMagnus3,kappas);
%             fprintf('\t Mean abs error Euler Ref and Magnus 3 at time %1.3g at center: %1.3e\n',...
%                             tk,errAbsEulerRefMagnus3{kk}(d,d));
            [errDisEulerRefMagnus3{kk},xEulerRefMagnus3{kk},fEulerRefMagnus3{kk},...
                nBlowupsEulerRefMagnus3{kk}]= ...
                errorDis(utEulerRef,utMagnus3,kappas,...
                        'blowup',blowup,...
                        'normType',errorDisType);
%             for iKappa=1:1:length(kappas)
%                 if nBlowupsEulerRefMagnus3{kk}{iKappa}>-1
%                     fprintf('\t\t Number of blowups for Magnus 3 with kappa=%d: %d\n',...
%                             kappas(iKappa),nBlowupsEulerRefMagnus3{kk}{iKappa})
%                 end
%             end
        end
        kk=kk+1;
    end
    fprintf('Total Elapsed time for Euler ref is %3.3g seconds.\n',ctimeEulerRefTotal)
    fprintf('Total Elapsed time for Magnus order 1 is %3.3g seconds.\n',ctimeMagnus1Total)
    fprintf('Total Elapsed time for Magnus order 2 is %3.3g seconds.\n',ctimeMagnus2Total)
    fprintf('Total Elapsed time for Magnus order 3 is %3.3g seconds.\n',ctimeMagnus3Total)
    %% Save results
    disp('Save results')
    fileName=sprintf('%s_T%1.3f_d%d_dtlog%1.2e_dt%1.3e_M%d_100%d%d%d',...
        methodName,T,nx,dtMagnusLog,dT,M,...
        useMagnus1,useMagnus2,useMagnus3);
    root=[pwd, '/' ,'Results'];
    matRoot=[root,'/','StepTest','/',methodName];
    mkDir(matRoot)
    if useMagnus1
        stepTest.Magnus1.ctime=ctimeMagnus1Total;
        stepTest.Magnus1=addAbsError(stepTest.Magnus1,'EulerRef',errAbsAvgEulerRefMagnus1{end});
        stepTest.Magnus1=addRelError(stepTest.Magnus1,'EulerRef',errDisEulerRefMagnus1{end});
        stepTest.Magnus1.T=T;
        stepTest.Magnus1.dt=dT;
        stepTest.Magnus1.nx=nx;
        stepTest.Magnus1.nv=nv;
        stepTest.Magnus1.M=M;
        stepTest.Magnus1.dtLog=dtMagnusLog;
        stepTest.Magnus1.kappas=kappas;
    end
    if useMagnus2
        stepTest.Magnus2.ctime=ctimeMagnus2Total;
        stepTest.Magnus2=addAbsError(stepTest.Magnus2,'EulerRef',errAbsAvgEulerRefMagnus2{end});
        stepTest.Magnus2=addRelError(stepTest.Magnus2,'EulerRef',errDisEulerRefMagnus2{end});
        stepTest.Magnus2.T=T;
        stepTest.Magnus2.dt=dT;
        stepTest.Magnus2.nx=nx;
        stepTest.Magnus2.nv=nv;
        stepTest.Magnus2.M=M;
        stepTest.Magnus2.dtLog=dtMagnusLog;
        stepTest.Magnus2.kappas=kappas;
    end
    if useMagnus3
        stepTest.Magnus3.ctime=ctimeMagnus3Total;
        stepTest.Magnus3=addAbsError(stepTest.Magnus3,'EulerRef',errAbsAvgEulerRefMagnus3{end});
        stepTest.Magnus3=addRelError(stepTest.Magnus3,'EulerRef',errDisEulerRefMagnus3{end});
        stepTest.Magnus3.T=T;
        stepTest.Magnus3.dt=dT;
        stepTest.Magnus3.nx=nx;
        stepTest.Magnus3.nv=nv;
        stepTest.Magnus3.M=M;
        stepTest.Magnus3.dtLog=dtMagnusLog;
        stepTest.Magnus3.kappas=kappas;
    end
    delFile([matRoot,'/',fileName,'.mat'])
    save([matRoot,'/',fileName,'.mat'],'stepTest')
end
disp('Done')
function sct=addAbsError(sct,ref,absError)
    sct.(ref).absError=absError;
end
function sct=addRelError(sct,ref,relError)
    temp=zeros(1,length(relError));
    for it1 = 1:1:length(temp)
        temp(it1)=mean(relError{it1},1);
    end
    sct.(ref).relError=temp;
end
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