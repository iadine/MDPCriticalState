function [LP, LPplus,LPminus] = mdp_critical_state(P, R, discount, precision)
% written by IC

% mdp_critical_state  Calculates Loss of performance when implementing a
% suboptimal action in a state
% Arguments -------------------------------------------------------------
% Let S = number of states, A = number of actions
%   P(SxSxA)  = transition matrix 
%               P could be an array with 3 dimensions 
%   R(SxSxA) or (SxA) = reward matrix
%              R could be an array with 3 dimensions (SxSxA) or 
%   discount  = discount rate in ]0; 1[
% Evaluation -------------------------------------------------------------
%   LPplus(S) = max loss value function, associated to a specific action
%   LPminus(S) = min loss value function, associated to a specific action
%--------------------------------------------------------------------------


% Perform value iteration
[sol_policy,~,~] = mdp_value_iteration(P, R, discount, precision);

% Evaluate policy iteratively
[y_Vpolicy,y_Q] = mdp_eval_policy_iterative_q(P, R, discount, sol_policy);

% Get the dimensions of Q
[nbS, nbA] = size(y_Q);

% Initialize LP, MinusLP, and PlusLP
LP = zeros(nbS, nbA);
LPminus = zeros(nbS, 1);
LPplus = zeros(nbS, 1);

% Calculate LP, MinusLP, and PlusLP
for i = 1:nbS
    % For all states
    optA = sol_policy(i);
    SetA = setdiff(1:nbA, optA);
    LP(i, optA) = NaN;
    for a = SetA
        LP(i, a) = (y_Vpolicy(i) - y_Q(i, a)) / y_Vpolicy(i);
    end
    LPminus(i) = min(LP(i, :), [], 'omitnan');
    LPplus(i) = max(LP(i, :), [], 'omitnan');
end



