function utMagnus=magnusConst(timegrid,W,A,B,phiVec,nx,nv,Ti,order,varargin)
%%MAGNUSCONST computes the Magnus expansion for the stochastic Langevin
% equation with constant coefficients.
%   Input:
%       A (nx*nv x nx*nv sparse array): Ito coefficient
%       B (nx*nv x nx*nv sparse array): Lebesgue coefficient
%       phiVec (nx*nv x M array): initial datum
%       nx (int): dimension of x grid
%       nv (int): dimension of v grid
%       Ti (1xNmagnus int array): evaluation points
%       varargin (cell array): name value pairs
%           'device' (default: gpu): 'cpu', 'gpu'
%           'precision' (default: half): 'double', 'single', 'half'
%   Output:
%
device='gpu';
precision = 'half';
for iV=1:1:length(varargin)
    switch varargin{iV}
        case 'device'
            device=varargin{iV+1};
        case 'precision'
            precision=varargin{iV+1};
    end
end

m=30;
N=length(timegrid);
M=size(W,4);
dt=timegrid(end)/(N-1);
timegrid=reshape(timegrid,[N, 1]);
W=reshape(W,[N,M]);

utMagnus=zeros(nx*nv,length(Ti)*M);
% switch device
%     case 'gpu'
%         LASTN = maxNumCompThreads(2);
% end
switch order
    case 1
        switch device
            case 'cpu'
                currT=timegrid(Ti,1);
                currW=W(Ti,:);
                parfor i=1:length(Ti)*M
                    [tiTemp,mi]=ind2sub([length(Ti),M],i);
                    O=firstorder(A,B,currT(tiTemp),currW(tiTemp,mi));
                   [F,~,~,~] = expmvtay2(O,phiVec(:,mi),m,precision);
                    utMagnus(:,i)=F;
                end
            case 'gpu'
                A=gpuArray(A);
                B=gpuArray(B);
                currT=timegrid(Ti,1);
                currW=W(Ti,:);
                parfor i=1:length(Ti)*M
                    [tiTemp,mi]=ind2sub([length(Ti),M],i);
                    O=firstorder(A,B,currT(tiTemp),currW(tiTemp,mi));
                    [F,~,~,~] = expmvtay2(O,gpuArray(phiVec(:,mi)),m,precision);
                    utMagnus(:,i)=gather(F);
                end
            case 'gpuLowMem'
                currT=timegrid(Ti,1);
                currW=W(Ti,:);
                parfor i=1:length(Ti)*M
                    [tiTemp,mi]=ind2sub([length(Ti),M],i);
                    O=firstorder(A,B,currT(tiTemp),currW(tiTemp,mi));
                    [F,~,~,~] = expmvtay2(gpuArray(O),gpuArray(phiVec(:,mi)),m,precision);
                    utMagnus(:,i)=gather(F);
                end
        end
    case 2
        IW=l_int(W,dt);
        switch device
            case 'cpu'
                BA=comm(B,A);
                currW=W(Ti,:);
                currT=timegrid(Ti,1);
                currIW=IW(Ti,:);
                parfor i=1:length(Ti)*M
                    [tiTemp,mi]=ind2sub([length(Ti),M],i);
                    O=firstorder(A,B,currT(tiTemp),currW(tiTemp,mi))+...
                        secondorder(A,BA,currT(tiTemp,1),currW(tiTemp,mi),currIW(tiTemp,mi));
                   [F,~,~,~] = expmvtay2(O,phiVec(:,mi),30,precision);
                    utMagnus(:,i)=F;
                end
            case 'gpu'
                A=gpuArray(A);
                B=gpuArray(B);
                currT=timegrid(Ti,1);
                currW=W(Ti,:);
                currIW=IW(Ti,:);
                BA=comm(B,A);
                parfor i=1:length(Ti)*M
                    [tiTemp,mi]=ind2sub([length(Ti),M],i);
                    O=firstorder(A,B,currT(tiTemp),currW(tiTemp,mi))+...
                        secondorder(A,BA,currT(tiTemp,1),currW(tiTemp,mi),currIW(tiTemp,mi));
                    [F,~,~,~] = expmvtay2(O,gpuArray(phiVec(:,mi)),m,precision);
                    utMagnus(:,i)=gather(F);
                end
            case 'gpuLowMem'
                currT=timegrid(Ti,1);
                currW=W(Ti,:);
                currIW=IW(Ti,:);
                BA=comm(B,A);
                parfor i=1:length(Ti)*M
                    [tiTemp,mi]=ind2sub([length(Ti),M],i);
                    O=firstorder(A,B,currT(tiTemp),currW(tiTemp,mi))+...
                        secondorder(A,BA,currT(tiTemp,1),currW(tiTemp,mi),currIW(tiTemp,mi));
                    [F,~,~,~] = expmvtay2(gpuArray(O),gpuArray(phiVec(:,mi)),m,precision);
                    utMagnus(:,i)=gather(F);
                end
        end
    case 3
        IW=l_int(W,dt);
        IsW=l_int(W.*timegrid,dt);
        IW2=l_int(W.^2,dt);
        switch device
            case 'cpu'
                BA=comm(B,A);
                BAA=comm(BA,A);
                BAB=comm(BA,B);
                currT=timegrid(Ti,1);
                currW=W(Ti,:);
                currIW=IW(Ti,:);
                currIsW=IsW(Ti,:);
                currIW2=IW2(Ti,:);
                parfor i=1:length(Ti)*M
                    [tiTemp,mi]=ind2sub([length(Ti),M],i);
                    O=firstorder(A,B,currT(tiTemp),currW(tiTemp,mi))+...
                        secondorder(A,BA,currT(tiTemp,1),currW(tiTemp,mi),currIW(tiTemp,mi))+...
                        thirdorder(BAA,BAB,currT(tiTemp,1),currW(tiTemp,mi),currIsW(tiTemp,mi),currIW(tiTemp,mi),currIW2(tiTemp,mi));
                    [F,~,~,~] = expmvtay2(O,phiVec(:,mi),m,precision);
                    utMagnus(:,i)=F;
                end
            case 'gpu'
                A=gpuArray(A);
                B=gpuArray(B);
                BA=comm(B,A);
                BAA=comm(BA,A);
                BAB=comm(BA,B);
%                 currT=gpuArray(timegrid(Ti,1));
%                 currW=gpuArray(W(Ti,:));
%                 currIW=gpuArray(IW(Ti,:));
%                 currIsW=gpuArray(IsW(Ti,:));
%                 currIW2=gpuArray(IW2(Ti,:));
                currT=timegrid(Ti,1);
                currW=W(Ti,:);
                currIW=IW(Ti,:);
                currIsW=IsW(Ti,:);
                currIW2=IW2(Ti,:);
                parfor i=1:length(Ti)*M
                    [tiTemp,mi]=ind2sub([length(Ti),M],i);
                    O=firstorder(A,B,currT(tiTemp),currW(tiTemp,mi))+...
                        secondorder(A,BA,currT(tiTemp,1),currW(tiTemp,mi),currIW(tiTemp,mi))+...
                        thirdorder(BAA,BAB,currT(tiTemp,1),currW(tiTemp,mi),currIsW(tiTemp,mi),currIW(tiTemp,mi),currIW2(tiTemp,mi));
                    [F,~,~,~] = expmvtay2(O,gpuArray(phiVec(:,mi)),m,precision);
                    utMagnus(:,i)=gather(F);
                end
            case 'gpuLowMem'
                BA=comm(B,A);
                BAA=comm(BA,A);
                BAB=comm(BA,B);
                currT=timegrid(Ti,1);
                currW=W(Ti,:);
                currIW=IW(Ti,:);
                currIsW=IsW(Ti,:);
                currIW2=IW2(Ti,:);
                parfor i=1:length(Ti)*M
                    [tiTemp,mi]=ind2sub([length(Ti),M],i);
                    O=firstorder(A,B,currT(tiTemp),currW(tiTemp,mi))+...
                        secondorder(A,BA,currT(tiTemp,1),currW(tiTemp,mi),currIW(tiTemp,mi))+...
                        thirdorder(BAA,BAB,currT(tiTemp,1),currW(tiTemp,mi),currIsW(tiTemp,mi),currIW(tiTemp,mi),currIW2(tiTemp,mi));
                    [F,~,~,~] = expmvtay2(gpuArray(O),gpuArray(phiVec(:,mi)),m,precision);
                    utMagnus(:,i)=gather(F);
                end
        end
%     case 4
%         switch device
%             case 'cpu'
% 
%             case 'gpu'
%         end
%         IW=l_int(W,dt);
%         IsW=l_int(W.*timegrid,dt);
%         IW2=l_int(W.^2,dt);
%         BA=comm(B,A);
%         BAA=comm(BA,A);
%         BAB=comm(BA,B);
%         BAAA=comm(BAA,A);
%         BABA=comm(BAB,A);
%         BABB=comm(BAB,B);
%         BAA2=comm(BA,A^2);
%         IsW2=l_int(W.^2.*timegrid,dt);
%         Is2W=l_int(W.*timegrid.^2,dt);
%         IW3=l_int(W.^3,dt);
%         parfor i=1:length(Ti)*M
%             [tiTemp,mi]=ind2sub([length(Ti),M],i);
%             ti=Ti(tiTemp);
%             O=firstorder(A,B,timegrid(ti,1),W(ti,mi))+...
%                 secondorder(A,BA,timegrid(ti,1),W(ti,mi),IW(ti,mi))+...
%                 thirdorder(BAA,BAB,timegrid(ti,1),W(ti,mi),IsW(ti,mi),IW(ti,mi),IW2(ti,mi))+...
%                 fourthorder(BAAA,BABA,BABB,BAA2,timegrid(ti,1),W(ti,mi),IsW(ti,mi),IsW2(ti,mi),Is2W(ti,mi),IW(ti,mi),IW2(ti,mi),IW3(ti,mi));
%             [F,~,~,~,~]=expmv(1,O,phiVec,[],precision);
%             utMagnus(:,i)=F;
%         end
    otherwise
        error('Order must be 1, 2, 3, 4: %d given.',order)
end
% switch device
%     case 'gpu'
%         LASTN = maxNumCompThreads(LASTN);
% end
utMagnus=reshape(utMagnus,nx,nv,length(Ti),M);
end
function O1=firstorder(A,B,timegrid,W)
    O1=B.*timegrid+...
        A.*W;
end
function O2=secondorder(A,BA,timegrid,W,IW)
    O2=-A*A.*.5.*timegrid+...
        BA.*IW-BA.*timegrid.*W.*.5;
end
function O3=thirdorder(BAA,BAB,timegrid,W,IsW,IW,IW2)
    O3=-BAB.*(1/12).*timegrid.^2.*W+...
      BAA.*(1/12).*timegrid.*...
        W.^2+...
      BAB.*...
        IsW-...
      BAA.*.5.*...
        W.*IW-...
      BAB.*.5.*timegrid.*...
        IW+...
      BAA.*.5.*...
        IW2;
end
% function O1=firstorder(A,B,timegrid,W)
%     O1=B.*timegrid+...
%         A.*W;
% end
% function O2=secondorder(A,BA,timegrid,W,IW)
%     O2=-A^2./2.*timegrid+...
%         BA.*IW-BA.*timegrid.*W./2;
% end
% function O3=thirdorder(BAA,BAB,timegrid,W,IsW,IW,IW2)
%     O3=-BAB./12.*timegrid.^2.*W+...
%       BAA./12.*timegrid.*...
%         W.^2+...
%       BAB.*...
%         IsW-...
%       BAA./2.*...
%         W.*IW-...
%       BAB./2.*timegrid.*...
%         IW+...
%       BAA./2.*...
%         IW2;
% end
% function O4=fourthorder(BAAA,BABA,BABB,BAA2,timegrid,W,IsW,IsW2,Is2W,IW,IW2,IW3)
%     O4=BAAA.*(IW.*W.^2./12-IW2.*W./4+IW3./6)+...
%         BABA.*(-IsW.*W./2+IsW2./2+IW.*timegrid.*W./6+W.^2.*timegrid.^2./24-IW2.*timegrid./4)+...
%         BABB.*(-IsW.*timegrid./2+Is2W./2+IW.*timegrid.^2./12)+...
%         BAA2.*(W.*timegrid.^2./24-IsW./2+IW.*timegrid./4);
% end
function I=l_int(W,dt)
    I=zeros(size(W));
    I(2:1:end,:)=cumsum(W(1:1:end-1,:).*dt,1); 
end
function Z=comm(A,B)
    Z=A*B-B*A;
end