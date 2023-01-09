function utExact=exactApprox(T,Ti,timegrid,xgrid,vgrid,a,sigma,W,phi)

    xigrid=linspace(-4,4,400)';
    etagrid=linspace(-4,4,400);
    dxi=(xigrid(end)-xigrid(1))/(length(xigrid)-1);
    deta=(etagrid(end)-etagrid(1))/(length(etagrid)-1);
    N=size(W,3);
    M=size(W,4);
    W=reshape(W,[N,M]);
    dt=T/(N-1);
    IW=l_int(W,dt);

%     Gamma0 = @(t,x,v) ...
%         (sqrt(3)/(pi .* t.^2 * (a-sigma.^2))).*...
%         exp(-2./(a-sigma.^2).*(v.^2./t - 2.*v.*x./t.^2 + 3.*x.^2./t.^3));
%     m1 = @(xi,eta,ti,wi) xi+timegrid(ti).*eta-sigma.*IW(ti,wi);
%     m2 = @(xi,eta,ti,wi) eta-sigma*W(ti,wi);
    phi = matlabFunction(phi,'Vars',[sym('x'),sym('v')]);
    Phi=phi(xigrid,etagrid);
    utExact=zeros(length(xgrid)*length(vgrid)*length(Ti)*M,1);
    currT=timegrid(Ti);
    currW=W(Ti,:);
    currIW=IW(Ti,:);
    parfor j=1:length(xgrid)*length(vgrid)*length(Ti)*M
        [xj,vj,tj,wj]=ind2sub([length(xgrid),length(vgrid),length(Ti),M],j);
        m1=xigrid+currT(tj).*etagrid-sigma.*currIW(tj,wj);
        m2=etagrid-sigma.*currW(tj,wj);
        currX=xgrid(xj)-m1;
        currV=vgrid(vj)-m2;
        Gamma=(sqrt(3)/(pi .* currT(tj).^2 * (a-sigma.^2))).*...
        exp(-2./(a-sigma.^2).*...
        (currV.^2./currT(tj) -...
         3.*currV.*currX./currT(tj).^2 +...
         3.*currX.^2./currT(tj).^3));
        utExact(j)=sum(Gamma.*Phi,[1,2]).*dxi.*deta;
    end
    utExact=reshape(utExact,length(xgrid),length(vgrid),length(Ti),M);
end
function I=l_int(W,dt)
    I=zeros(size(W));
    I(2:1:end,:)=cumsum(W(1:1:end-1,:).*dt,1); 
end