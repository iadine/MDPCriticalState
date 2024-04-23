% function state= seeState(index,state_matrix)
% Action:
%   return state at index 'index'
% input 
%   index = index for getting its state
% output 
%   state = state at index {state for Abalone,state for sea otter}
% author:   iadine.chades@csiro.au
%
function state= seeState(index,state_matrix)

state=zeros(1,2);
% get state
state(1)=state_matrix(index,2);
state(2)=state_matrix(index,3);
