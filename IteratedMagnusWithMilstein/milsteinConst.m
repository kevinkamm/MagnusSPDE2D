function ut=milsteinConst(TTi,T,Neuler,dWeuler,Dx,Dv,Dxx,Dvv,H,Fx,Fv,Gxx,Gxv,Gvv,Sx,Sv,phi,varargin)
%%EULER computes the Euler-Maruyama scheme for the stochastic Langevin
% equation with constant coefficients.
% Magnus (A, B).
%   Input:
%       D (nv x nv sparse array): -v part in SPDE
%       E (nx x nx sparse array): first derivative wrt x (central diff)
%       F (nv x nv sparse array): first derivative wrt v (central diff)
%       G (nv x nv sparse array): second derivative wrt v (central diff)
%       varargin (cell array): name value pairs
%           'Device' (default: cpu): 'cpu', 'gpu'
%           'Mode' (default: full): 'sparse', 'full'
%   Output:

device='cpu';
mode='full';
if isempty(varargin)==0
    for k=1:2:length(varargin)
        switch varargin{k}
            case 'device'
                device=varargin{k+1};
            case 'mode'
                mode=varargin{k+1};
        end
    end
end
M=size(dWeuler,4);
nx=size(Dx,1);
nv=size(Dv,1);
dt=T/(Neuler-1);
switch mode
    case 'full'
        Dx=full(Dx);
        Dv=full(Dv);
        Dxx=full(Dxx);
        Dvv=full(Dvv);
        H=full(H);
        Fx=full(Fx);
        Fv=full(Fv);
        Gxx=full(Gxx);
        Gxv=full(Gxv);
        Gvv=full(Gvv);
        Sx=full(Sx);
        Sv=full(Sv);
    case 'sparse'
    otherwise
        error("Unsupported mode: " +...
               mode + " given; full or sparse expected");
end
switch device
    case 'gpu'
        Dx=gpuArray(Dx);
        Dv=gpuArray(Dv);
        Dxx=gpuArray(Dxx);
        Dvv=gpuArray(Dvv);
        H=gpuArray(H);
        Fx=gpuArray(Fx);
        Fv=gpuArray(Fv);
        Gxx=gpuArray(Gxx);
        Gxv=gpuArray(Gxv);
        Gvv=gpuArray(Gvv);
        Sx=gpuArray(Sx);
        Sv=gpuArray(Sv);
%         dWeuler=gpuArray(dWeuler);
%         ut=gpuArray.zeros(nx,nv,Neuler,M);
%         ut(:,:,1,:)=repmat(gpuArray(phi),[1,1,1,M]);
        utTemp=zeros(nx,nv,1,M);
        ut=zeros(nx,nv,length(TTi),M);
        utTemp(:,:,1,:)=phi;
    case 'cpu'
        utTemp=zeros(nx,nv,1,M);
        ut=zeros(nx,nv,length(TTi),M);
        utTemp(:,:,1,:)=phi;
    otherwise
        error("Unsupported device: " +...
               device + " given; cpu or gpu expected");
end
% O=ones(nx,nv);
% Dsigma=(Sx.*leftMult(Dx,O)+Sv.*rightMult(O,Dv'));
% Dsigma=(Sx.*Dx+Sv.*Dv');
ti=1;
for tk=1:1:Neuler-1
    if tk==TTi(ti)
        ut(:,:,ti,:)=utTemp;
        ti=ti+1;
    end
    sigma=(Sx.*leftMult(Dx,utTemp)+Sv.*rightMult(utTemp,Dv'));
    Dsigma=(Sx.*leftMult(Dx,sigma)+Sv.*rightMult(sigma,Dv'));
    utTemp=...
        utTemp+...
        (H.*utTemp+...
        Fx.*leftMult(Dx,utTemp)+...
        Fv.*rightMult(utTemp,Dv')+...
        (Gxx./2).*leftMult(Dxx,utTemp)+...
        Gxv.*rightMult(leftMult(Dx,utTemp),Dv')+...
        (Gvv./2).*rightMult(utTemp,Dvv')).*dt+...
        sigma.*dWeuler(:,:,tk,:)+...
        Dsigma.*(dWeuler(:,:,tk,:).^2-dt)./2;
end
if Neuler==TTi(ti)
    ut(:,:,ti,:)=utTemp;
end
function Z=leftMult(A,X)
    switch string(device)+" "+string(mode)
        case "gpu full"
            X=gpuArray(X);
%             Z=pagefun(@mtimes,A,X); 
            Z=pagemtimes(A,X);
        case "cpu full"
            Z=pagemtimes(A,X);
        case "cpu sparse"
            Z=zeros(nx,nv,1,M);
            for mi=1:1:M
                Z(:,:,1,mi)=A*X(:,:,mi);
            end
        case "gpu sparse"
            X=gpuArray(X);
            Z=gpuArray.zeros(nx,nv,1,M);
            for mi=1:1:M
                Z(:,:,1,mi)=A*X(:,:,mi);
            end
        otherwise
            error('Unknown device or mode')
    end
end
function Z=rightMult(X,A)
    switch string(device)+" "+string(mode)
        case "gpu full"
            X=gpuArray(X);
%             Z=pagefun(@mtimes,X,A); 
            Z=pagemtimes(X,A);
        case "cpu full"
            Z=pagemtimes(X,A);
        case "cpu sparse"
            Z=zeros(nx,nv,1,M);
            for mi=1:1:M
                Z(:,:,1,mi)=X(:,:,mi)*A;
            end
        case "gpu sparse"
            X=gpuArray(X);
            Z=gpuArray.zeros(nx,nv,1,M);
            for mi=1:1:M
                Z(:,:,1,mi)=X(:,:,mi)*A;
            end
        otherwise
            error('Unknown device or mode')
    end
end
end