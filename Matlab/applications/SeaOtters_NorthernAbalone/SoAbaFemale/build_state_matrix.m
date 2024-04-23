% function state_matrix= build_state_matrix(nbs_aba,nbs_seaotter)
% Action: 
%    Builds the index state matrix that provides the corresponding state
%    number for a given abalone state and sea otter state; and vice versa
% Input:
%    nbs_aba: max number of abalone states
%    nbs_seaotter: max number of sea otter states
% Output:
%   state_matrix = matrix (index, state aba, state so)
% Side effect:
%   can be quite large depending on the discretisation used.
% Author: iadine.chades@csiro.au

function state_matrix= build_state_matrix(nbs_aba,nbs_seaotter)

nbs=nbs_aba*nbs_seaotter;
state_matrix=zeros(nbs,3);
state_matrix(:,1)=[1:nbs];

% corresponding state
for k=1:nbs_seaotter
    state_matrix(((k-1)*nbs_aba+1):k*nbs_aba,1:3)=...
        [((k-1)*nbs_aba+1):k*nbs_aba;0:(nbs_aba-1);ones(1,nbs_aba)*(k-1)]';
end;
