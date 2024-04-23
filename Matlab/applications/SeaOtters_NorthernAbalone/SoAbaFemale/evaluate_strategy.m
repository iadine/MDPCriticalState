% function evaluate_strategy(policyPi,str,iFR,iFRnum)
% Action:
%   simulate the dynamic of sea otters and abalone with a stochastic
%   occurence of oil spills (*) assuming a functional response and a
%   management strategy is applied (policyPi). This function also provides
%   figures 3-4 of the paper Chades et al.(in review)
% Input:
%   policyPi = a management strategy that maps an action for each index
%   state
%   str = strings that identifies strategy and FR
%   iFR = FR family
%   iFRnum = FR number
% Output:
%   kickass figure saved in directory Fig/
%
% Side effect:
%   Best to close figures as you go to avoid memory issue.
%
% Author: iadine.chades@csiro.au

function evaluate_strategy(policyPi,str,iFR,iFRnum)

global PARAM_QL PARAM_ABALONE
global IS_DISPLAYED_GRAPH

if IS_DISPLAYED_GRAPH==1

% Initializing parameters
T=PARAM_QL.Time_Horizon;    % Optimisation time
t=PARAM_QL.t;               % Decision time step
sea_o=zeros(T,t);       
m=zeros(T,t);
t_aba=zeros(T,t);
decision=zeros(T,1);
Oil=zeros(T,t);
TAba=zeros(T*t+1,1);
Aba100mm=zeros(T*t+1,1);

figure('color','white','name',['Example of a simulation for ',str]);
SoPop=0;
[AbaPop,AbaPopF]=initialising_northern_abalone(0);    % 0 = no plot
init_state=[abaloneDensity2state(sum(AbaPop)/PARAM_ABALONE.area),...
            SOabundance2state(SoPop)];  % Corresponding MDP state
current_state=init_state;           % internal variable
TAba(1)=sum(AbaPop);                % Tracking the whole abalone population
Aba100mm(1,:)=sum(AbaPop(6:10));    % Tracking the adults >100mm    

for i=1:T
    action=policyPi(seeIndex(current_state));
    if (action==1 && current_state(2) ~=0)
        disp('GRRR cant introduce more SO');
    end;
    if action==-1
        disp('0_o conservation action not well defined!');
        disp(action);
    end
    decision(i)=action;

    [AbaPop,AbaPopF,TAbaAdults,TSoPop,next_state,outcome,oilSpill]=simulation_t(AbaPop,AbaPopF,SoPop,action,iFR,iFRnum);
    m(i,:)=outcome';
    sea_o(i,:)=TSoPop';
    Aba100mm(i+1,:)=sum(AbaPop(6:10));
    SoPop=TSoPop(end);
    t_aba(i,:)=TAbaAdults';
    Oil(i,:)=oilSpill';
    current_state=next_state;
end

OS=zeros(T*t+1,1);
So=zeros(T*t+1,1);
Rm=zeros(T*t+1,1);

% Trick for graphs useful when t<>1
for k=2:T+1
    OS(k+(k-2)*(t-1):k+(k-2)*(t-1)+(t-1))=Oil(k-1,1:t)';
    So(k+(k-2)*(t-1):k+(k-2)*(t-1)+(t-1))=sea_o(k-1,1:t)';
    TAba(k+(k-2)*(t-1):k+(k-2)*(t-1)+(t-1))=t_aba(k-1,1:t)';
    Rm(k+(k-2)*(t-1):k+(k-2)*(t-1)+(t-1))=m(k-1,1:t)';    
end

subplot(4,1,1:2);
% Plotting aba pop, whole population and 100mm population
x=0:T*t;

h(1)=line(x,TAba/PARAM_ABALONE.area,'color',[1 0.5 0],'LineStyle','-');
hold on
h(4)=line(x,Aba100mm/PARAM_ABALONE.area,'color','r','LineStyle','-');
hold on

% Trick to get SO pop on a different axis
ax1=gca;
set(ax1,'XColor','w','YColor',[1 0.5 0],'YLim',[0 1.5]);
set(get(ax1,'yLabel'),'String','Abalone density (m^{-2})');
ax2 = axes('Position',get(ax1,'Position'),...
    'XAxisLocation','bottom',...
    'YAxisLocation','right',...
    'Color','none',...
    'XColor','k','YColor',[0 0.5 1]);
set(get(ax2,'yLabel'),'String','Sea otter abundance');
hold on
h(2)=line(x,So,'color',[0 0.5 1],'LineStyle','--');
xlabel('Time (years)')
hold on

% Trick to plot oil spill
dec=0;
for k=2:T*t+1
    if OS(k)~=0
        if dec==0
            dec=100;
        else
            dec=0;
        end
        hold on
        h(3)=plot(ax2,k-1,So(k)+200, 'k*');
        hold off
    end;
end
hold on

% Ploting rewards
subplot(4,1,3);
sth=stem(x,Rm,'-k');
set(sth,'Marker','.')


hold on;
legend('Reward','location','best');
legend('boxoff');
box off

% Ploting actions
subplot(4,1,4);
A=[
    0.8 0.8 0.8;
    0.8 0 0.2;
    0 1 0.2;
    0.1 0.1 0.1
    0.5 0.5 0;
    ];
image(decision','CDataMapping','scaled')
axis off
colormap(A);
caxis([0 5]);
end_name='ONE_RUN';
DirFig='Fig/';
saveas(gcf,[DirFig,str,end_name], 'fig');

end
