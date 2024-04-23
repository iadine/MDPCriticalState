% function expi=compute_expert_policy(expert_num)
% Action:
%   Generates the expert strategies following:
%       expert_num = 1 {introduction of SO/ antipoaching}
%       expert_num = 2 {introduction of SO/ antipoaching/culling}
%       expert_num = 3 {introduction of SO/ antipoaching/ culling & antipoaching}
%       expert_num = 4 {introduction of SO/ do nothing}
%       if expert_num~{1..4} ERROR MESSAGE
% Input:
%   expert_num = [1,4]
% Output:
%   expi = the strategy that maps an action to an index state number.
% Side effects:
%   The difference between each strategy appears above the sea otter 
% parameter kculling threshold defined in load_parameter()
%
% Author: iadine.chades@csiro.au
%

function expi=compute_expert_policy(expert_num)
global PARAM_MDP PARAM_SO

limit_state_culling=SOabundance2state(PARAM_SO.kculling*PARAM_SO.k)+1; %first state we can do culling

nbs=PARAM_MDP.nbs_aba*PARAM_MDP.nbs_so;   % nbs = number of states
expi=zeros(nbs,1);

if expert_num==1
    for i=1:nbs
        s=seeState(i,PARAM_MDP.state_matrix);
        if s(2)==0        %no otter
            expi(i)=1;   %introduction of SO
        else
                expi(i)=2;  %Antipoach aba
        end
    end
elseif expert_num==2
    for i=1:nbs
        s=seeState(i,PARAM_MDP.state_matrix);
        if s(2)==0        %no otter
            expi(i)=1;   %introduction of SO
        else
            if s(2)>limit_state_culling+1
                expi(i)=3;  %culling SO
            else
                expi(i)=2;  %Antipoach aba
            end
        end
    end
elseif expert_num==3
    for i=1:nbs
        s=seeState(i,PARAM_MDP.state_matrix);
        if s(2)==0        %no otter
            expi(i)=1;   %introduction of SO
        else
            if s(2)>limit_state_culling+1
                expi(i)=4;  %culling SO+antipoach
            else
                expi(i)=2;  %Antipoach aba
            end
        end
    end
elseif expert_num==4
    for i=1:nbs
        s=seeState(i,PARAM_MDP.state_matrix);
        if s(2)==0        %no otter
            expi(i)=1;   %introduction of SO
        else
            if s(2)>limit_state_culling+1
                expi(i)=0;  %ne rien faire
            else
                expi(i)=0;  %ne rien faire
            end
        end
    end
else
    disp('Incorrect expert strategy number in compute_expert_policy');
    error('Expert strategy should be an integer between [1,4]');
end