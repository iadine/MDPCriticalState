% function main_SDP(nsim,                % number of simulation
%                   isAPoachEfficient,   % 0 or 1    reduce poaching by 0.5/0.75 (default 0)
%                   isDisplayOn,         % 0 or 1    figures off/on   (default 0)
%                   maxKAba,             % carrying capacity in best habitat (default 3.34)
%                   maxrAba,             % max growth rate in best habitat (default 1.6)
%                   KSO,                 % carrying capacity sea otter (default 4073)
%                   rSO)                 % growth rate sea otter (default 0.191)
%
% main_SDP is the main function to run the optimisation procedure.
%
% Action:
% This files computes the SDP solutions to the sea otter and abalone
% problem. Can SO and NA co-exist? by Chades et al (under review)
%
% Input: 
%   nsim,                % number of simulation
%   isAPoachEfficient,   % 0 or 1    reduce poaching by 0.5/0.75 (default 0)
%   isDisplayOn,         % 0 or 1    figures off/on   (default 0)
%   maxKAba,             % carrying capacity in best habitat (default 3.34)
%   maxrAba,             % max growth rate in best habitat (default 1.6)
%   KSO,                 % carrying capacity sea otter (default 4073)
%   rSO)                 % growth rate sea otter (default 0.191)
% Output:
%   Run this file to compute the SDP solutions, evaluate the performances by
% simulations.
% - The simulation results are stored in Results/[date]_results_VI.txt
% - The figures are saved in Fig/
%
% Note that this procedure takes a long time ~12hours-24 hours as.
%
% DISCLAIMER: this is a research program, not suitable for commercial
% purpose. Please contact the author for updates and fixed bugs.
%
% AUTHOR: iadine chades iadine.chades@csiro.au
% Last update: June 2011
% Matlab R2009b

function main_SDP(varargin)
rand('state', 300);     % Same random seed

%------- checking input arguments
numvarargs = length(varargin);
if numvarargs > 7
    error('myfuns:main_SDP:TooManyInputs', ...
        'requires at most 7 optional inputs');
end
% set defaults for optional inputs
optargs = {50 0 0 3.34 1.6 4073 0.191};
optargs(1:numvarargs) = varargin;
% Place optional args in memorable variable names
[nsim, isAPoachEfficient, isDisplayOn, maxKAba, maxrAba, KSO, rSO] = optargs{:};
% ------- end checking input arguments.

disp('-|| Optimisation procedure:')
disp('-|| Can sea otters and northern abalone co-exist under best management?')
disp('-|| Chades et al (under review)')

global Tr
global gfilename
global DIR_results
global IS_DISPLAYED_GRAPH
global PARAM_LINEAR_FR

IS_DISPLAYED_GRAPH=isDisplayOn;   % 1= yes; 0=No

% Initialising and loading parameters
disp('-> Loading parameters');
load_param(isAPoachEfficient,maxKAba,maxrAba,KSO,rSO);
disp('<- Parameters loaded');

% To reduce computation time use a small value for ksim (<100)
ksim=nsim;     % Number of simulation to estimate Tr
k_max=500;   % Number of simulation to evaluate the SDP strategy

% File management
output_filename=[DIR_results,date,'_results_VI']; % For example: 02-Jun-2011_results_
gfilename=fopen([output_filename,'.txt'],'w+');
fprintf(gfilename,'%d %f %f %f %f\n',isAPoachEfficient,maxKAba,...
        maxrAba,KSO,rSO);  % saving some common parameters with the file
fprintf(gfilename,'VI \n| Average rew | +- | SO abundance | +- | Aba density | +- | Aba100mm | +- ');
% End file management

nFR=3;  % There are 3 families of functional responses (FR):
% sigmoid, hyperbolic and linear.

for iFR=1:nFR
    if iFR==1       % Sigmoid FR
        beginTag='Sig_';
        for iSig=1:2     % 2 sigmoid FR
            debutt=[beginTag,num2str(iSig)];
            disp(['<START> Optimisation procedure for functional response ',debutt]);
            Tr=compute_transition(iFR,iSig,ksim);   % Estimation of Tr
            Policy_VI=value_iteration();            % Optimisation algorithm
            set_policy=explore_Q(Policy_VI,0.99);   % Similarities
            if IS_DISPLAYED_GRAPH==1, draw_policy(Policy_VI,...
                    [debutt,'_VI'],set_policy); end
            if IS_DISPLAYED_GRAPH==1, evaluate_strategy(Policy_VI,...
                    [debutt,'VI'],iFR,iSig); end
            evaluate_strategy_mean(Policy_VI,[debutt,'VI'],iFR,iSig,k_max);
            disp(['<\END> Optimisation procedure for ',debutt]);
        end
    elseif iFR==2   % Hyperbolic FR
        beginTag='Hyp_';
        for iHyp=1:2    % 2 hyperbolic FR
            debutt=[beginTag,num2str(iHyp)];
            disp(['<START> Optimisation procedure for functional response ',debutt]);
            Tr=compute_transition(iFR,iHyp,ksim);     % Estimation of Tr
            Policy_VI=value_iteration();       % Optimisation algorithm
            set_policy=explore_Q(Policy_VI,0.99); % Similarities
            if IS_DISPLAYED_GRAPH==1, draw_policy(Policy_VI,...
                    [debutt,Pi_best_filename],set_policy); end
            if IS_DISPLAYED_GRAPH==1, evaluate_strategy(Policy_VI,...
                    [debutt,'VI'],iFR,iHyp); end
            evaluate_strategy_mean(Policy_VI,[debutt,'VI'],iFR,iHyp,k_max);
            disp(['<\END> Optimisation procedure for ',debutt]);
        end;
    elseif iFR==3   % Linear FR
        beginTag='Linear_';
        for iLinear=1:length(PARAM_LINEAR_FR) % variable number of Linear FR
            debutt=[beginTag,num2str(iLinear)];
            disp(['<START> Optimisation procedure for functional response ',debutt]);
            Tr= compute_transition(iFR,iLinear,ksim);   % Estimation of Tr
            Policy_VI= value_iteration();               % Optimisation algorithm
            set_policy= explore_Q(Policy_VI,0.99);      % Similarities
            if IS_DISPLAYED_GRAPH==1, draw_policy(Policy_VI,...
                    [debutt,Pi_best_filename],set_policy); end
            if IS_DISPLAYED_GRAPH==1, evaluate_strategy(Policy_VI,...
                    [debutt,'VI'],iFR,iLinear); end
            evaluate_strategy_mean(Policy_VI,[debutt,'VI'],iFR,iLinear,k_max);
            disp(['<\END> Optimisation procedure for ',debutt]);
        end
    end
end
fclose(gfilename);
disp(['Performance results saved in ', output_filename,'.txt'])
disp('-|| END Optimisation procedure')

end