% function rewardSum=evaluate_strategy_mean(PolicyPi,str,iFR,iFRnum,nsim)
% Action:
%   evaluation by simulation of the dynamic of sea otters and abalone 
%   assuming a functional response and a management strategy is applied 
%   (policyPi). The average performances are saved as a line in
%   'date'-results files. This function also provides supplementary figures 
%   of the paper Chades et al.(in review).
% Input:
%   policyPi = a management strategy that maps an action for each index
%   state
%   str = strings that identifies strategy and FR
%   iFR = FR family
%   iFRnum = FR number
%   nsim   = number of simulations (>1000)
% Output:
%   kickass figure and result saved in directories Fig/ and Results/
% Global: gfilename
% Side effect:
%   Best to close figures as you go to avoid memory issue.
%
% Author: iadine.chades@csiro.au
%
function rewardSum=evaluate_strategy_mean(PolicyPi,str,iFR,iFRnum,nsim)

global PARAM_QL PARAM_ABALONE G 
global gfilename IS_DISPLAYED_GRAPH

T=PARAM_QL.Time_Horizon;
t=PARAM_QL.t;

klim=nsim;
rewardSum=zeros(klim,1);

Am=zeros(klim,T*PARAM_QL.t);
%Astd=zeros(klim,PARAM_QL.t);

SOm=zeros(klim,T*PARAM_QL.t);
%SOstd=zeros(klim,PARAM_QL.t);

Aba100mm=zeros(klim,T*PARAM_QL.t+1,1);
SoPop=0;
[AbaPop,AbaPopF]=initialising_northern_abalone(0); % 0 = no plot
init_state=[abaloneDensity2state(sum(AbaPop)/PARAM_ABALONE.area),...
            SOabundance2state(SoPop)];  % Corresponding MDP state

index = seeIndex(init_state);
Init=struct('AbaPop',AbaPop,'SoPop',SoPop,'G',G,'Init_state',init_state,'index',index);

for k=1:klim
    % Mean of klim simulations
    rSum=0;
    SoPop=0;
    in_aba=Init.AbaPop;
    AbaPop=Init.AbaPop;
    init_state=Init.Init_state;
    current_state=init_state;
    Aba100mm(k,1)=sum(AbaPop(6:10));
    for L=1:T
        actionS=PolicyPi(seeIndex(current_state));
        if actionS==-1
            display('Incorrect action number')
        end;
        [AbaPop,AbaPopF,TAbaAdults,TSoPop,next_state,outcome]=...
            simulation_t(AbaPop,AbaPopF,SoPop,actionS,iFR,iFRnum);
        Aba100mm(k,L+1)=sum(AbaPop(6:10)); % aba plus de 100 mmm
        SoPop=TSoPop(end);
        rSum=rSum+sum(outcome);
        current_state=next_state;

        Am(k,L*t-t+1:L*t)=TAbaAdults;
        SOm(k,L*t-t+1:L*t)=TSoPop;
    end
    rewardSum(k)=rSum;
end

ci = [0.025 0.975]; % Confidence levels
cii = floor(ci*(klim-1))+1;

sr = sort(Am,1);
NewAm=mean(sr(cii(1):cii(2),:));
NewAstd=std(sr(cii(1):cii(2),:));
NTAba=[sum(in_aba), NewAm];
NAstd=[0, NewAstd];

% average line 100mm Abalone
sr = sort(Aba100mm,1);
NewA100=mean(sr(cii(1):cii(2),:));
NewA100std=std(sr(cii(1):cii(2),:));
NTAba100=NewA100;
NAstd100=NewA100std;


sr = sort(SOm,1);
NewSO=mean(sr(cii(1):cii(2),:));
NewSOstd=std(sr(cii(1):cii(2),:));
NSo=[0,NewSO];
NSostd=[0,NewSOstd];


%% display results
if IS_DISPLAYED_GRAPH ==1
    figure('color','white','name',['Average simulation results for ',str]);
    subplot(4,2,1:4);
    x=0:T*t;
    xx=0:T*t;

    % CI
    hold on
    nh1=line(x,NTAba/PARAM_ABALONE.area,'color',[1 0.5 0],'LineStyle','-');
    hold on
    nh11=line(x,(NTAba+NAstd)/PARAM_ABALONE.area,'color',[1 0.5 0],'Marker','.');
    hold on
    nh12=line(x,(NTAba-NAstd)/PARAM_ABALONE.area,'color',[1 0.5 0],'Marker','.');
    nth1=text(5,1.3,num2str(NTAba(T*t+1)/PARAM_ABALONE.area,3));

    hold on
    h100mm=line(xx,NTAba100/PARAM_ABALONE.area,'color',[1 0. 0.1],'LineStyle','-');
    nth1=text(5,1.1,num2str(NTAba100(T+1)/PARAM_ABALONE.area,3));
    % end CI

    hold on
    ax1=gca;
    set(ax1,'XColor','w','YColor',[1 0.5 0],'YLim',[0 1.4],'visible','on');
    set(get(ax1,'yLabel'),'String','Abalone density (m^{-2})');
    ax2 = axes('Position',get(ax1,'Position'),...
        'XAxisLocation','bottom',...
        'YAxisLocation','right',...
        'Color','none',...
        'XColor','k','YColor',[0 0.5 1]);
    set(get(ax2,'yLabel'),'String','Sea otter abundance');
    hold on

    %CI
    nh2=line(x,NSo,'color',[0 0.5 1],'LineStyle','-');
    hold on
    nh21=line(x,NSo+NSostd,'color',[0 0.5 1],'Marker','.');
    hold on
    nh22=line(x,NSo-NSostd,'color',[0 0.5 1],'Marker','.');
    xlabel('Time (years)')
    hold on
    %nth2=text(T*t-20,NSo(T*t+1),num2str(NSo(T*t+1),3));
    % end CI

    %set(nth2,'visible','on')
    set(nth1,'visible','on')
    set(nh1,'visible','on')
    set(nh11,'visible','on')
    set(nh12,'visible','on')


    subplot(4,2,5:8);
    show_FR(iFR,iFRnum);
    
    end_name='AVERAGE';
    DirFig='Fig/';
    saveas(gcf,[DirFig,str,end_name], 'fig');
end

% SAVE RESULTS in FILE gfilename
fprintf(gfilename,'\n%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f', ...
    mean(rewardSum),std(rewardSum),NSo(T*t),...
    NSostd(T*t),NTAba(T*t)/PARAM_ABALONE.area,NAstd(T*t)/PARAM_ABALONE.area,...
    NTAba100(T*t+1)/PARAM_ABALONE.area,NAstd100(T*t+1)/PARAM_ABALONE.area );
end