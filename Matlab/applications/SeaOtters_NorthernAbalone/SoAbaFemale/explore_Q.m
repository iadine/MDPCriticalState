% function liste_a=explore_Q(PiB,sensitivity)
% Action: 
%   explore how far the actions are from each other given a Q matrix
%   computed from an optimisation procedure. This function is meant to
%   provide some information in case other actions might be closed to be
%   "near-optimal" (based on the sensitivity parameter).
% Input:
%   PiB: a strategy computed
%   sensitivity: number between 0.9 and 1.
% Output: 
%    liste_a: list of actions so that q(s,list_a(s))/q(s,PiB(s))>
%    sensitivity
% Global:
%   q: is defined as a global variable
% Side effect:
%   This function is particularly interesting to understand the strategies
%   computed using reinforcement learning procedure, where some actions 
%   might have very similar output.
%
% Author: iadine.chades@csiro.au

function liste_a=explore_Q(PiB,sensitivity)


global q PARAM_MDP
global limit_state_culling

nbs=PARAM_MDP.nbs_aba*PARAM_MDP.nbs_so;   % nbs = number of states
liste_a=zeros(nbs,4)*-1;

for i=1:nbs
    s=seeState(i,PARAM_MDP.state_matrix);
    if s(2)==0
        TAction=[1,2];
    elseif s(2)<limit_state_culling
        TAction=[2];
    elseif s(2)>=limit_state_culling
        TAction=[2,3,4];
    end
    for a=1:length(TAction)
            if q(i,PiB(i))~=0 && q(i,TAction(a))/q(i,PiB(i))>sensitivity
                liste_a(i,TAction(a))=q(i,TAction(a))/q(i,PiB(i));
            end
    end
end
