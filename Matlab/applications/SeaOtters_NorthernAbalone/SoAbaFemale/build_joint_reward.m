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
function R=build_joint_reward(k_aba,level_so,level_aba)
global  IS_DISPLAYED_GRAPH

x=0:0.005:k_aba;
y=interp1(level_aba,level_so,x,'pchip');
gap=(100)/(size(x,2)-1);
z=0:gap:100;
R=[x',y',z'];

if IS_DISPLAYED_GRAPH==1
figure('color','white');
plot(level_aba,level_so,'o',x,y)
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

end
