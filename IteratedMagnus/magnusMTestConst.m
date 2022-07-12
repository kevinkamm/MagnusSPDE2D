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
useExact = true;
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
[h,fx,fv,gxx,gxv,gvv,sx,sv,phi,hasExact,methodName]=...
    selectCoefficients('exact1',[a,sigma]);
% [h,fx,fv,gxx,gxv,gvv,sx,sv,phi,hasExact,methodName]=...
%     selectCoefficients('Langevin1',[a,sigma]);
%% Video & Plot parameters
backgroundColor = 'w';
% backgroundColor = [53,54,58]./255;
textColor = 'k';
% textColor = [237,237,237]./255;

%% Simulations
% M = 100;
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
dT=0.05;
dtMagnusLog = 1e-4;
M=[10:10:90,100:100:1000];
for M = flip(M)
    %% Initialize variables
    ctimeMagnus1Total=0;
    ctimeMagnus2Total=0;
    ctimeMagnus3Total=0;
    ctimeExactTotal=0;
    errDisExactMagnus1={};
    errDisExactMagnus2={};
    errDisExactMagnus3={};
    errAbsExactMagnus1={};
    errAbsExactMagnus2={};
    errAbsExactMagnus3={};
    errAbsAvgExactMagnus1={};
    errAbsAvgExactMagnus2={};
    errAbsAvgExactMagnus3={};

    Nmagnus = floor(dT/dtMagnusLog)+1;
    phiMagnus1=repmat(Phi(:),1,M);
    phiMagnus2=repmat(Phi(:),1,M);
    phiMagnus3=repmat(Phi(:),1,M);
    
    kk=1;
    WmagnusCell=cell(floor(T/dT),1);
    WexactCell=cell(floor(T/dT),1);
%     disp('Generate Brownian motion')
    for tk=dT:dT:T
        [~,~,Wmagnus,~,~,~,~]=brownianMotion(dT,Nmagnus,[],[],M);
        WmagnusCell{kk}=Wmagnus;
        WexactCell{kk}=Wmagnus;
        kk=kk+1;
    end
    kk=1;
    fprintf('Curr dT=%g\n',dT)
    for tk=dT:dT:T
        %% Brownian motion
        tnMagnus=Nmagnus;
        Wmagnus=WmagnusCell{kk};
        
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
        %% Calculate Exact
        tExact=linspace(0,dT,Nmagnus);
%         disp('Compute Exact')
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
%         fprintf('Elapsed time for Exact is %3.3g seconds.\n',ctimeExact)
        %% Errors
%         disp('Reference Method exact')
        if useMagnus1
            [errAbsExactMagnus1{kk},errAbsAvgExactMagnus1{kk}]=...
                absError(utExact,utMagnus1,kappas);
%             fprintf('\t Mean abs error Exact and Magnus 1 at time %1.3g at center: %1.3e\n',...
%                         tk,errAbsExactMagnus1{kk}(d,d));
            [errDisExactMagnus1{kk},xExactMagnus1{kk},fExactMagnus1{kk},...
                nBlowupsExactMagnus1{kk}]= ...
                errorDis(utExact,utMagnus1,kappas,...
                        'blowup',blowup,...
                        'normType',errorDisType);
%             for iKappa=1:1:length(kappas)
%                 if nBlowupsExactMagnus1{kk}{iKappa}>-1
%                     fprintf('\t\t Number of blowups for Magnus 1 with kappa=%d: %d\n',...
%                             kappas(iKappa),nBlowupsExactMagnus1{kk}{iKappa})
%                 end
%             end
        end
        if useMagnus2
            [errAbsExactMagnus2{kk},errAbsAvgExactMagnus2{kk}]=...
                absError(utExact,utMagnus2,kappas);
%             fprintf('\t Mean abs error Exact and Magnus 2 at time %1.3g at center: %1.3e\n',...
%                         tk,errAbsExactMagnus2{kk}(d,d));
            [errDisExactMagnus2{kk},xExactMagnus2{kk},fExactMagnus2{kk},...
                nBlowupsExactMagnus2{kk}]= ...
                errorDis(utExact,utMagnus2,kappas,...
                        'blowup',blowup,...
                        'normType',errorDisType);
%             for iKappa=1:1:length(kappas)
%                 if nBlowupsExactMagnus2{kk}{iKappa}>-1
%                     fprintf('\t\t Number of blowups for Magnus 2 with kappa=%d: %d\n',...
%                             kappas(iKappa),nBlowupsExactMagnus2{kk}{iKappa})
%                 end
%             end
        end
        if useMagnus3
            [errAbsExactMagnus3{kk},errAbsAvgExactMagnus3{kk}]=...
                absError(utExact,utMagnus3,kappas);
%             fprintf('\t Mean abs error Exact and Magnus 3 at time %1.3g at center: %1.3e\n',...
%                         tk,errAbsExactMagnus3{kk}(d,d));
            [errDisExactMagnus3{kk},xExactMagnus3{kk},fExactMagnus3{kk},...
                nBlowupsExactMagnus3{kk}]= ...
                errorDis(utExact,utMagnus3,kappas,...
                        'blowup',blowup,...
                        'normType',errorDisType);
%             for iKappa=1:1:length(kappas)
%                 if nBlowupsExactMagnus3{kk}{iKappa}>-1
%                     fprintf('\t\t Number of blowups for Magnus 3 with kappa=%d: %d\n',...
%                             kappas(iKappa),nBlowupsExactMagnus3{kk}{iKappa})
%                 end
%             end
        end
        kk=kk+1;
    end
    fprintf('Total Elapsed time for Magnus order 1 is %3.3g seconds.\n',ctimeMagnus1Total)
    fprintf('Total Elapsed time for Magnus order 2 is %3.3g seconds.\n',ctimeMagnus2Total)
    fprintf('Total Elapsed time for Magnus order 3 is %3.3g seconds.\n',ctimeMagnus3Total)
    fprintf('Total Elapsed time for Exact is %3.3g seconds.\n',ctimeExactTotal)
    %% Save results
    disp('Save results')
    fileName=sprintf('%s_T%1.3f_d%d_dtlog%1.2e_dt%1.3e_M%d_100%d%d%d',...
        methodName,T,nx,dtMagnusLog,dT,M,...
        useMagnus1,useMagnus2,useMagnus3);
    root=[pwd, '/' ,'Results'];
    matRoot=[root,'/','MTest'];
    mkDir(matRoot)
    if useMagnus1
        mTest.Magnus1.ctime=ctimeMagnus1Total;
        mTest.Magnus1=addAbsError(mTest.Magnus1,'Exact',errAbsAvgExactMagnus1{end});
        mTest.Magnus1=addRelError(mTest.Magnus1,'Exact',errDisExactMagnus1{end});
        mTest.Magnus1.T=T;
        mTest.Magnus1.dt=dT;
        mTest.Magnus1.nx=nx;
        mTest.Magnus1.nv=nv;
        mTest.Magnus1.M=M;
        mTest.Magnus1.dtLog=dtMagnusLog;
        mTest.Magnus1.kappas=kappas;
    end
    if useMagnus2
        mTest.Magnus2.ctime=ctimeMagnus2Total;
        mTest.Magnus2=addAbsError(mTest.Magnus2,'Exact',errAbsAvgExactMagnus2{end});
        mTest.Magnus2=addRelError(mTest.Magnus2,'Exact',errDisExactMagnus2{end});
        mTest.Magnus2.T=T;
        mTest.Magnus2.dt=dT;
        mTest.Magnus2.nx=nx;
        mTest.Magnus2.nv=nv;
        mTest.Magnus2.M=M;
        mTest.Magnus2.dtLog=dtMagnusLog;
        mTest.Magnus2.kappas=kappas;
    end
    if useMagnus3
        mTest.Magnus3.ctime=ctimeMagnus3Total;
        mTest.Magnus3=addAbsError(mTest.Magnus3,'Exact',errAbsAvgExactMagnus3{end});
        mTest.Magnus3=addRelError(mTest.Magnus3,'Exact',errDisExactMagnus3{end});
        mTest.Magnus3.T=T;
        mTest.Magnus3.dt=dT;
        mTest.Magnus3.nx=nx;
        mTest.Magnus3.nv=nv;
        mTest.Magnus3.M=M;
        mTest.Magnus3.dtLog=dtMagnusLog;
        mTest.Magnus3.kappas=kappas;
    end
    delFile([matRoot,'/',fileName,'.mat'])
    save([matRoot,'/',fileName,'.mat'],'mTest')
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