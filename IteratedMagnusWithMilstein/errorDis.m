function [errorDis,x,f,varargout]=errorDis(utRef,utApprox,kappas,varargin)
    [nx,nv,M] = size(utRef,[1,2,4]);
    normType='rel1'; 
    percentile=100000;
    blowup=0;
    if ~isempty(varargin)
        for k=1:2:length(varargin)
            switch varargin{k}
                case 'normType'
                    normType=varargin{k+1};
                case 'percentile'
                    percentile=varargin{k+1};
                case 'blowup'
                    blowup=varargin{k+1};
            end
        end
    end
    errorDis=cell(length(kappas),1);
    x=cell(length(kappas),1);
    f=cell(length(kappas),1);
    nBlowups=cell(length(kappas),1);
    for iKappa = 1:1:length(kappas)
        kappa=2^iKappa; 
        regionX=1+floor(nx/2-nx/(2*kappa)):1:ceil(nx/2+nx/(2*kappa));
        regionV=1+floor(nv/2-nv/(2*kappa)):1:ceil(nv/2+nv/(2*kappa));
        switch normType
            case 'rel1'
                tempDis=...
                    rel1(utRef(regionX,regionV,end,:),...
                          utApprox(regionX,regionV,end,:));
                if blowup
                    indNoBlowup=tempDis<1;
                    nBlowups{iKappa}=M-sum(indNoBlowup,'all');
                    tempDis=tempDis(indNoBlowup);
                end
                errorDis{iKappa}=tempDis;
            case 'rel2'
                tempDis=...
                    rel2(utRef(regionX,regionV,end,:),...
                          utApprox(regionX,regionV,end,:));
                if blowup
                    indNoBlowup=tempDis<1;
                    nBlowups{iKappa}=M-sum(indNoBlowup,'all');
                    tempDis=tempDis(indNoBlowup);
                end
                errorDis{iKappa}=tempDis;
            case 'abs1'
                tempDis=...
                    abs1(utRef(regionX,regionV,end,:),...
                          utApprox(regionX,regionV,end,:));
                if blowup
                    indNoBlowup=tempDis<1e6;
                    nBlowups{iKappa}=M-sum(indNoBlowup,'all');
                    tempDis=tempDis(indNoBlowup);
                end
                errorDis{iKappa}=tempDis;
        end
        %Empirical CDF at time t
        if ~isempty(tempDis)
            [tempF,tempX]=ecdf(tempDis);
            c=find(tempF>percentile,1,'first');
            if ~isempty(c)
                tempF=tempF(1:1:c-1);
                tempX=tempX(1:1:c-1);
            end
        else
            tempX=0;
            tempF=0;
        end
        x{iKappa}=tempX;
        f{iKappa}=tempF;
    end
    if nargout>3
        varargout{1}=nBlowups;
    end
    function Err_t=rel1(uRef,uApprox)
        Err_t=zeros(M,1);
        nom = squeeze(sqrt(sum(abs(uRef-uApprox).^2,[1,2])));
        denom = squeeze(sqrt(sum(uRef.^2,[1,2])));
        indNonZero=denom > 0;
        Err_t(indNonZero)=nom(indNonZero)./denom(indNonZero);
    end
    function Err_t=rel2(uRef,uApprox)
        Err_t=zeros(M,1);
        denom = uRef;
        nom = uApprox;
        indNonZero=denom > 0;
        Err_t(indNonZero)=mean(nom(indNonZero)./denom(indNonZero) -1,[1,2]);
    end
    function Err_t=abs1(uRef,uApprox)
        Err_t=mean(abs(uRef-uApprox),[1,2]);
    end
end