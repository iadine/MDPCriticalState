% Function seaotterState = SOabundance2state(popSeaotter)
% Action:
%   Convert abundance into a MDP state
% Input:
%   popSeaOtter = number of sea otters
% Output:
%   seaotterState = corresponding state
% Author: iadine.chades@csiro.au

function seaotterState = SOabundance2state(popSeaotter)
global PARAM_SO PARAM_MDP

K=PARAM_SO.k;
nbs=PARAM_MDP.nbs_so;    %21
tranche=K/(nbs-1);
if popSeaotter==0
    seaotterState=0;
else
    seaotterState=fix(popSeaotter/tranche)+1;
    if seaotterState>=nbs
        seaotterState=nbs-1;
    end
end

