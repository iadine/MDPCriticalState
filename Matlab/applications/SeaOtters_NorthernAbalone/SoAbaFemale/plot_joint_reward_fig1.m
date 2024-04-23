% function R=build_joint_reward(k_aba,level_so,level_aba)
% Action: 
%    Builds/prints the joint reward of the optimisation model. Given
%    thresholds for both species. Because we are interested in managing
%    both species at the same time, we must account for the species with
%    the lowest "status". This function will be used in the optimisation 
%    procedure so that the immediate reward for both species corresponds to 
%    R(r,aba)= min(r(so),r(aba)).
% Input:
%   level_so, level_aba: Recovery stages for each species
%   k_aba: carrying capacity abalone
% Global variable: 
%   IS_DISPLAYED_GRAPH = turns graph on/off
% Output:
%   R: a reward matrix (Abalone ; sea otters ; reward)
% Side effects: 
%   Unknown
% Author: iadine.chades@csiro.au
%
function R=plot_joint_reward_fig1()
global  IS_DISPLAYED_GRAPH
IS_DISPLAYED_GRAPH=1;
k_aba=1.92;
level_so= [0 0.3  0.4 0.6 1 ];
level_aba=[0 0.02 0.1 0.5 k_aba];

x=0:0.005:k_aba;
y=interp1(level_aba,level_so,x,'pchip');
gap=(100)/(size(x,2)-1);
z=0:gap:100;
R=[x',y',z'];


% points of interest corresponding ~ to the recovery stages
pts=[R(1,:);R(5,:);R(21,:);R(101,:);R(385,:);]

if IS_DISPLAYED_GRAPH==1
figure('color','white');
plot(level_aba,level_so,'x',x,y)
hold on
xlabel('Abalone density');
ylabel('%K Sea Otter');
title('Interpolated Joint Reward Function');
hold on; 
plot3(x,y,z,'--g')
zlabel('reward')
hold on
grid
view(-14,33);
box off
figure('color','white');
h(1)=line(x,z,'Color','r');
hold on
plot(pts(:,1),pts(:,3),'kx','MarkerSize',10,'LineWidth',2)
ylabel('Reward');
ax1=gca;
set(ax1,'Xlim',[0 k_aba]);
set(ax1,'XColor','r','YColor','k');
xlabel('Abalone density (per m^2)','Color','k');
ax2 = axes('Position',get(ax1,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','b','YColor','k');
hold on
xlabel('Sea otter abundance (% carrying capacity)','Color','k');
ylabel('Reward');
%ylabel('%K Sea Otter');
hold on
h(2)=line(y,z,'Color','b', 'Parent', ax2);
lh=legend(h,'Abalone reward function (JR_{aba})','Sea otter reward function (JR_{so})');
set(lh,'Box','off','Location','NorthWest')
plot(pts(:,2),pts(:,3),'kx','MarkerSize',10,'LineWidth',2)

end
