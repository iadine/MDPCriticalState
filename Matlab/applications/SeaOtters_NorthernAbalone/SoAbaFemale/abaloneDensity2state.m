% function abaloneState=abaloneDensity2state(d_aba)
% Action: 
%   Map an abalone state to a density 
% Input:
%   d_aba = density of abalone;
% Output:
%   abaloneState = state id (0 to nbstates-1)
% Side effects:
%   The last state might be larger than others due to d_aba going above k.
%
% Author: iadine.chades@csiro.au
%
function abaloneState=abaloneDensity2state(d_aba)

global PARAM_ABALONE PARAM_MDP

if d_aba<PARAM_ABALONE.k
    abaloneState=fix(d_aba/PARAM_ABALONE.discUnit); % aba state starts at 0 
elseif d_aba >= PARAM_ABALONE.k
    disp('density higher than or equal to k_aba');
    abaloneState=PARAM_MDP.nbs_aba -1;
end;
