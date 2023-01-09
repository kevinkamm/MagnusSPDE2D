function [A,B,Dx,Dv,Dxx,Dvv,H,Fx,Fv,Gxx,Gxv,Gvv,Sx,Sv,Phi]=coefficients(xgrid,vgrid,h,fx,fv,gxx,gxv,gvv,sx,sv,phi)
%%COEFFICIENTS computes the coefficients for Euler and Magnus.
%   Input:
%       xgrid (nx+2 x 1 array): position grid
%       vgrid (1 x nv+2 array): velocity grid
%       h (function handle, args x,v): drift coefficient for zero
%                                      derivative
%       fx (function handle, args x,v): drift coefficient for first 
%                                       derivative wrt x
%       fv (function handle, args x,v): drift coefficient for first 
%                                       derivative wrt v
%       gxx (function handle, args x,v): drift coefficient for second 
%                                        derivative wrt x
%       gvv (function handle, args x,v): drift coefficient for second 
%                                        derivative wrt x
%       gxv (function handle, args x,v): drift coefficient for second 
%                                        derivative wrt x and v
%       sx (function handle, args x,v): diffusion coefficient for first 
%                                       derivative wrt x
%       sv (function handle, args x,v): diffusion coefficient for first 
%                                       derivative wrt v
%   Output:
%       A (nx*nv x nx*nv sparse array): Magnus coefficient for Ito integral
%       B (nx*nv x nx*nv sparse array): Magnus coefficient for Lebesgue integral
%       Dx (nx x nx sparse array): first derivative wrt x (central diff)
%       Dv (nv x nv sparse array): first derivative wrt v (central diff)
%       Dxx (nv x nv sparse array): second derivative wrt x (central diff)
%       Dvv (nv x nv sparse array): second derivative wrt v (central diff)
%       H
%       Fx
%       Fv
%       Gxx
%       Gxv
%       Gvv
%       Sx
%       Sv
%
%   Usage:
%       [A,B,Dx,Dv,Dxx,Dvv,H,Fx,Fv,Gxx,Gxv,Gvv,Sx,Sv]=
%           coefficients(xgrid,vgrid,h,fx,fv,gxx,gxv,gvv,sx,sv)
%

nx=length(xgrid)-2; %subtract boundary conditions
nv=length(vgrid)-2; %subtract boundary conditions
dx=(xgrid(end)-xgrid(1))/(length(xgrid)-1);
dv=(vgrid(end)-vgrid(1))/(length(vgrid)-1);

Dx=0.5.*(1./dx).*spdiags([-1,0,1].*ones(nx,1),-1:1,nx,nx);
Dv=0.5.*(1./dv).*spdiags([-1,0,1].*ones(nv,1),-1:1,nv,nv);
Dxx=(1./(dx.^2)).*spdiags([1,-2,1].*ones(nx,1),-1:1,nx,nx);
Dvv=(1./(dv.^2)).*spdiags([1,-2,1].*ones(nv,1),-1:1,nv,nv);
H=broadcast(nx,nv,h(xgrid(2:end-1),vgrid(2:end-1)));
Fx=broadcast(nx,nv,fx(xgrid(2:end-1),vgrid(2:end-1)));
Fv=broadcast(nx,nv,fv(xgrid(2:end-1),vgrid(2:end-1)));
Gxx=broadcast(nx,nv,gxx(xgrid(2:end-1),vgrid(2:end-1)));
Gxv=broadcast(nx,nv,gxv(xgrid(2:end-1),vgrid(2:end-1)));
Gvv=broadcast(nx,nv,gvv(xgrid(2:end-1),vgrid(2:end-1)));
Sx=broadcast(nx,nv,sx(xgrid(2:end-1),vgrid(2:end-1)));
Sv=broadcast(nx,nv,sv(xgrid(2:end-1),vgrid(2:end-1)));
Phi=broadcast(nx,nv,phi(xgrid(2:end-1),vgrid(2:end-1)));
Ix=speye(nx,nx);
Iv=speye(nv,nv);
B=largeFullSparseMult(Fx(:),kron(Iv,Dx))+...
    largeFullSparseMult(Fv(:),kron(Dv,Ix))+...
    largeFullSparseMult(Gxx(:)./2,kron(Iv,Dxx))+...
    largeFullSparseMult(Gxv(:),kron(Dv,Dx))+...
    largeFullSparseMult(Gvv(:)./2,kron(Dvv,Ix));
if numel(H)==1
    H=H.*ones(nx,nv);
end
B=B+spdiags(H(:),1,nx*nv,nx*nv);
    
A=largeFullSparseMult(Sx(:),kron(Iv,Dx))+...
    largeFullSparseMult(Sv(:),kron(Dv,Ix));
end
function Y=broadcast(nx,nv,X)
%     if size(X,1)==1 && size(X,2)==1
%         Y=X.*ones(nx,nv);
    if size(X,1)==1 && size(X,2)==nv
        Y=ones(nx,1).*X;
    elseif size(X,1)==nx && size(X,2)==1
        Y=X.*ones(1,nv);
    else
        Y=X;
    end
end
function Y=largeFullSparseMult(v,X)
    if numel(v)==1
        Y=v.*X;
    else
%         Y=gather(distributed(v).*X);
        Y=spdiags(v,0,size(X,1),size(X,2))*X;
    end
end
