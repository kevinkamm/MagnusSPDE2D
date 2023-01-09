function [dWeulerRef,dWeuler,WeulerRef,Wmagnus,tiERM,tiERE,tiEM]=brownianMotion(T,NeulerRef,Neuler,Nmagnus,M)
%%BROWNIANMOTION computes the Brownian motion with different time steps for
% Euler and Magnus
%   Input:
%       T (1 x 1 double): the finite time horizon
%       NeulerRef (1 x 1 int): the number of time-steps for reference Euler
%       Neuler (1 x 1 int): the number of time-steps for Euler
%       Nmagnus (1 x 1 int): the number of time-steps for Magnus
%       M (1 x 1 int): the number of simulations
%   Output:
%       dWeuler (1 x 1 x Neuler-1 x M): increment of BM for Euler
%       tiEuler (1 x Nmagnus): indices to match Euler and Magnus time grids


    if ~isempty(Nmagnus)
        tiERM=1:1:Nmagnus;
        tiERM(2:1:end)=tiERM(1:1:end-1).*floor((NeulerRef-1)/(Nmagnus-1))+1;
    else
        tiERM=[];
    end

    if ~isempty(Neuler)
        tiERE=1:1:Neuler;
        tiERE(2:1:end)=tiERE(1:1:end-1).*floor((NeulerRef-1)/(Neuler-1))+1;
    else
        tiERE=[];
    end

    if ~isempty(Nmagnus) && ~isempty(Neuler)
        tiEM=1:1:Nmagnus;
        tiEM(2:1:end)=tiEM(1:1:end-1).*floor((Neuler-1)/(Nmagnus-1))+1;
    else
        tiEM=[];
    end
    
    N=randn(1,1,NeulerRef-1,M);
    
    dtEulerRef=T/(NeulerRef-1);
    dWeulerRef=sqrt(dtEulerRef).*N;
    
    WeulerRef=zeros(1,1,NeulerRef,M);
    WeulerRef(1,1,2:end,:)=dWeulerRef;
    WeulerRef=cumsum(WeulerRef,3);

    if ~isempty(Neuler)
        dWeuler=diff(WeulerRef(1,1,tiERE,:),1,3);
    else
        dWeuler=[];
    end
%     dWmagnus=diff(WeulerRef(1,1,tiERM,:),1,3);
    if ~isempty(Nmagnus)
        Wmagnus=WeulerRef(1,1,tiERM,:);
    else
        Wmagnus=[];
    end
%     %% test indices
%     tEulerRef=linspace(0,T,NeulerRef);
%     tEuler=linspace(0,T,Neuler);
%     tMagnus=linspace(0,T,Nmagnus);
%     
%     sum(abs(tEulerRef(tiERE)-tEuler),'all')
%     sum(abs(tEulerRef(tiERM)-tMagnus),'all')
%     
%     %% test paths
%     WeulerRef=zeros(1,1,NeulerRef,M);
%     WeulerRef(1,1,2:end,:)=dWeulerRef;
%     WeulerRef=cumsum(WeulerRef,3);
%     
%     Weuler=zeros(1,1,Neuler,M);
%     Weuler(1,1,2:end,:)=dWeuler;
%     Weuler=cumsum(Weuler,3);
%     
%     Wmagnus=zeros(1,1,Nmagnus,M);
%     Wmagnus(1,1,2:end,:)=dWmagnus;
%     Wmagnus=cumsum(Wmagnus,3);
%     
%     sum(abs(WeulerRef(1,1,tiERE,:)-Weuler),'all')
%     sum(abs(WeulerRef(1,1,tiERM,:)-Wmagnus),'all')
%     
%     %% test BM by plot
%     figure();hold on;
%     w=1;
%     plot(tEulerRef,squeeze(WeulerRef(1,1,:,w)),'r-');
%     plot(tEuler,squeeze(Weuler(1,1,:,w)),'bo');
%     plot(tMagnus,squeeze(Wmagnus(1,1,:,w)),'gx');
end

% function [dWeuler,Weuler,Wmagnus,tiEuler,tiMagnus]=brownianMotion(T,Neuler,Nmagnus,M)
% %%BROWNIANMOTION computes the Brownian motion with different time steps for
% % Euler and Magnus
% %   Input:
% %       T (1 x 1 double): the finite time horizon
% %       Neuler (1 x 1 int): the number of time-steps for Euler
% %       Nmagnus (1 x 1 int): the number of time-steps for Magnus
% %       M (1 x 1 int): the number of simulations
% %   Output:
% %       dWeuler (1 x 1 x Neuler-1 x M): increment of BM for Euler
% %       Wmagnus (1 x 1 x Nmagnus x M): BM for Magnus
% %       tiEuler (1 x Nmagnus): indices to match Euler and Magnus time grids
%     
%     Nmax=max(Neuler,Nmagnus);
%     dtMax=T/(Nmax-1);
%     dWmax=zeros(1,1,Nmax,M);
%     dWmax(1,1,2:end,:)=sqrt(dtMax).*randn(Nmax-1,M);
%     Wmax=cumsum(dWmax,3);
%     if Nmax==Neuler
%         tiMagnus=1:1:Nmagnus;
%         tiEuler=1:1:Nmagnus;
%         tiEuler(2:1:end)=tiEuler(1:1:end-1).*floor((Neuler-1)/(Nmagnus-1))+1;
%         Weuler=Wmax(1,1,:,:);
%         dWeuler=dWmax(1,1,2:end,:);
%         Wmagnus=Wmax(1,1,tiEuler,:);
%     else
%         tiEuler=1:1:Neuler;
%         tiMagnus=1:1:Neuler;
%         tiMagnus(2:1:end)=tiMagnus(1:1:end-1).*floor((Nmagnus-1)/(Neuler-1))+1;
%         Wmagnus=Wmax(1,1,:,:);
%         dWeuler=diff(Wmax(1,1,tiMagnus,:),1,3);
%         Weuler=Wmax(1,1,tiMagnus,:);
%     end
%     
% %     %% test
% %     Wtest=zeros(1,1,Nmagnus,M);
% %     Wtest(1,1,2:end,:)=diff(Wmax(1,1,tiEuler,:),1,3);
% % %     Wtest(1,1,2:end,:)=Wmax(1,1,tiMagnus(2:end),:)-Wmax(1,1,tiMagnus(1:end-1),:);
% %     Wtest=cumsum(Wtest,3);
% %     sum(abs(Wtest - Wmagnus),'all')
% end