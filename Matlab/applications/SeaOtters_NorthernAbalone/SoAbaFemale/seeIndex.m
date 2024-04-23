
% function index= seeIndex(state)
% Action:
%   mapping function
% Input: 
%    state = abalone and sea otter states {Abalone, Sea otter}
% Output 
%   index= index # of state 'state'
% Author:
%   iadine.chades@csiro.au

function index= seeIndex(state)

global PARAM_MDP
index=(PARAM_MDP.nbs_aba)*state(2)+state(1)+1;
end


