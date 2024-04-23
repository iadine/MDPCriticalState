% script abalone_study_adults() - Generate figure S3 of Chades et al (2012)
% Example: abalone_study_adults()
% Source: Chad�s I., Curtis J.M.R, Martin T.G. (2012) Setting realistic 
% recovery targets for interacting endangered species. Under review, Conservation Biology. 
% Author: iadine chades iadine.chades@csiro.au
% Copyright: See README.txt

function abalone_study_adults_female()
% Global variables

global PARAM_ABALONE PARAM_QL PARAM_POACHING

% Initialising/loading parameters
disp('-> Loading parameters');
load_param(0,3.34,1.6,4073,0.191);
disp('<- Parameters loaded');

figure('color','white');

nb_simu=1000; % change if needed

%% stochastic with NO poaching
T=PARAM_QL.Time_Horizon;
AbaPop=zeros(nb_simu,T);
NAbaPop=zeros(nb_simu,T);
for j=1:nb_simu
    [~,N1f]=initialising_northern_abalone(0,PARAM_ABALONE.init_pop);

    for i=1:T
        [N1,N1f]=northern_abalone_growth_t(N1f);
        AbaPop(j,i)=sum(N1f(6:end))/PARAM_ABALONE.area;
    end
end

xlabel('Time (years)')
ylabel('Adult abalone density per m��(>100mm)');

% sort the pop and remove 2.5% of the best and worse trajectories
[~,IX]=sort(AbaPop(:,T));
nb_r=round(nb_simu*0.025);
for i=1:nb_simu 
    NAbaPop(i,:)=AbaPop(IX(i),:);
end
NAbaPop=NAbaPop(nb_r:nb_simu-nb_r,:);
s = std(NAbaPop(:,:),1,1);
M=mean(NAbaPop);

line(1:T,M,'color','k','LineStyle','-');
text(T,M(T),num2str(M(T),'%2.2f'));
hold on
line(1:T,M+s,'color','k','LineStyle','-.');
hold on
line(1:T,M-s,'color','k','LineStyle','-.');
hold on
disp('>>Poaching Low =')
disp(PARAM_POACHING.low) 
disp('++Average density of adult ABALONE and +-SD 95%CI =')
disp(M(T))
disp(s(T))


%% stochastic with low poaching
T=PARAM_QL.Time_Horizon;
AbaPop=zeros(nb_simu,T);
NAbaPop=zeros(nb_simu,T);
for j=1:nb_simu
    [~,N1f]=initialising_northern_abalone(0,PARAM_ABALONE.init_pop);
    for i=1:T
        [N1,N1f]=northern_abalone_growth_t(N1f);
        [N1,N1f]=compute_poaching_impact(N1,N1f,0);% low poaching
        AbaPop(j,i)=sum(N1f(6:end))/PARAM_ABALONE.area;
    end
end

% order the pop and remove 2.5% of the best and worse trajectories
[~,IX]=sort(AbaPop(:,T));
nb_r=round(nb_simu*0.025);
for i=1:nb_simu
    NAbaPop(i,:)=AbaPop(IX(i),:);
end
NAbaPop=NAbaPop(nb_r:nb_simu-nb_r,:);
s = std(NAbaPop(:,:),1,1);
M=mean(NAbaPop);

line(1:T,M,'color','r','LineStyle','-');
text(T,M(T),num2str(M(T),'%2.2f'));
hold on
line(1:T,M+s,'color','r','LineStyle','-.');
hold on
line(1:T,M-s,'color','r','LineStyle','-.');
hold on
disp('>>Poaching Low =')
disp(PARAM_POACHING.low) 
disp('++Average density of adult ABALONE and +- SD 95%CI =')
disp(M(T))
disp(s(T))
%% stochastic with high med poaching

Aba=zeros(nb_simu,T);
NAba=zeros(nb_simu,T);
for j=1:nb_simu
    [~,N1f]=initialising_northern_abalone(0,PARAM_ABALONE.init_pop);
    for i=1:T
        [N1,N1f]=northern_abalone_growth_t(N1f);
        [N1,N1f]=compute_poaching_impact(N1,N1f,1);% high to medium poaching
        Aba(j,i)=sum(N1f(6:end))/PARAM_ABALONE.area;
    end
end

% order the pop and remove 2.5% of the best and worse trajectories
[r,IX]=sort(Aba(:,T));
for i=1:nb_simu
    NAba(i,:)=Aba(IX(i),:);
end
NAba=NAba(nb_r:nb_simu-nb_r,:);

s = std(NAba(:,:),1,1);
M=mean(NAba);
line(1:T,M,'color','k','LineStyle','-');
text(T,M(T),num2str(M(T),'%2.2f'));
hold on
line(1:T,M+s,'color','k','LineStyle','-.');
%text(T,M(T)+s(T),num2str(M(T)+s(T)));
hold on
line(1:T,M-s,'color','k','LineStyle','-.');
%text(T,M(T)-s(T),num2str(M(T)-s(T)));
%ylim([0 PARAM_ABALONE.k]);

disp('**!!Poaching high-med')
disp(PARAM_POACHING.high) 
disp(PARAM_POACHING.med)
disp('** Average density of adult ABALONE and +-SD 95%CI')
disp(M(T))
disp(s(T))
