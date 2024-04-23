function [Xmean,ci,A]=sea_otters_study()
%% Sea otter population model including risk of oil spill
% Xmean : population moyenne avec oil spill
% A : population sans oil spill
% l'aspect stochastique vient de l'occurence de marée noire.

global PARAM_SO PARAM_QL

% Initialising/loading parameters
disp('-> Loading parameters');
load_param(0,3.34,1.6,4073,0.191);
disp('<- Parameters loaded');

T=PARAM_QL.Time_Horizon;
nb_simu=1000;
K=PARAM_SO.k;

X=zeros(nb_simu,T);
A=zeros(1,T);
NewX=zeros(nb_simu,T);
OS=zeros(nb_simu,T);
for i=1:nb_simu
    N=PARAM_SO.init_pop;
    for j=1:T
        [X(i,j),OS(i,j)]=sea_otter_growth(N);
        N=X(i,j);
    end
end

nb_r=round(nb_simu*0.025); % nombre a enlever
[r,IX]=sort(X(:,T));
for i=1:nb_simu
    NewX(i,:)=X(IX(i),:);
end

NewX=NewX(nb_r:nb_simu-nb_r,:); % j'enleve au debut et a la fin
Xmean=mean(NewX);
Xstd=std(NewX);
ci=Xstd(T);

%% sans maree noire

PARAM_SO.os_freq=0; %no oil spill
N=PARAM_SO.init_pop;
for j=1:T
    [A(j),B]=sea_otter_growth(N);
    N=A(j);
end

figure('color','w');
% plot 
line(1:T,Xmean,'color','k','linestyle','-');
hold on;
line(1:T,Xmean+Xstd,'color','k','linestyle','.');
hold on;
line(1:T,Xmean-Xstd,'color','k','linestyle','.');
hold on;
line(1:T,A,'color','g','linestyle','-');

% line (0:T,ones(1,T+1)*0.3*K,'color',[0.6 0.6 0.6],'linestyle','--');
% text(3,0.25*K,'Endangered','HorizontalAlignment','left');
% hold on
% line (0:T,ones(1,T+1)*0.4*K,'color',[0.6 0.6 0.6],'linestyle','--');
% text(3,0.35*K,'Threatened','HorizontalAlignment','left');
% hold on
% line (0:T,ones(1,T+1)*0.6*K,'color',[0.6 0.6 0.6],'linestyle','--');
% text(3,0.55*K,'Special concern','HorizontalAlignment','left');
% text(3,0.75*K,'Not a Risk','HorizontalAlignment','left');
% 
text(T+1,Xmean(T),num2str(Xmean(T),'%5.1f\n'));
text(T+1,Xmean(T)+Xstd(T),num2str(Xmean(T)+Xstd(T),'%5.1f\n'));
text(T+1,Xmean(T)-Xstd(T),num2str(Xmean(T)-Xstd(T),'%5.1f\n'));

% ---> change that
% less than 60% of K => Status = Endangered 60%-80% of K => Status =
% Threatened 80% K => Status = sensitive from (Abalone ? Urchin ? Sea Otter
% Meeting) plot states transition line

xlabel('Time (year)');
ylabel('Sea otter abundance');
%title(['Sea otter population model including oil spill (frequency=',num2str(oil_spill_frequency),') - K=',num2str(K) ,' r=',num2str(r)])
% legend(gca,'BVH','BVH+std','BVH-std','no oil spill','location','southOutside')
% legend(gca,'boxoff')
box off

PARAM_SO.os_freq=0.1; % on remet a jour la variable globale