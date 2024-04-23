% function main_experts(nsim,               % simulation number
%                       isAPoachEfficient,   % 0 or 1    reduce poaching by 0.5/0.75 (default 0)
%                       isDisplayOn,         % 0 or 1    figures off/on   (default 0)
%                       maxKAba,             % carrying capacity in best habitat (default 3.34)
%                       maxrAba,             % max growth rate in best habitat (default 1.6)
%                       KSO,                 % carrying capacity sea otter (default 4073)
%                       rSO)                 % growth rate sea otter (default 0.191))
%
% Action: main_experts is the main function to simulate the performances of
% the predifined strategies. main_expert simulates and evaluates the 
% predefined strategies (4) under different functional responses 
% (linear, sigmoid, hyperbolic)
%
% Input: 
%   nsim,                % number of simulations (default 500)
%   isAPoachEfficient,   % 0 or 1    reduce poaching by 0.5/0.75 (default 0)
%   isDisplayOn,         % 0 or 1    figures off/on   (default 0)
%   maxKAba,             % carrying capacity in best habitat (default 3.34)
%   maxrAba,             % max growth rate in best habitat (default 1.6)
%   KSO,                 % carrying capacity sea otter (default 4073)
%   rSO)                 % growth rate sea otter (default 0.191)%   
%   other parameters are initialized in function load_param.
%
% Output:  
% i) The simulation results are stored in
% Results/[date]_results_expert[1-4].txt
% ii) The figures are saved in Fig/ 
% 
% You can use results_analysis to print out results for paper.
%
% DISCLAIMER: this is a research program, not suitable for commercial
% purpose. Please contact the author for updates.
% 
% AUTHOR: iadine chades iadine.chades@csiro.au
% Last update: June 2011
% Matlab R2009b

function main_experts(varargin)

%------- checking input arguments
numvarargs = length(varargin);
if numvarargs > 7
    error('myfuns:main_experts:TooManyInputs', ...
        'requires at most 6 optional inputs');
end
% set defaults for optional inputs
optargs = {500 0 0 3.34 1.6 4073 0.191};
optargs(1:numvarargs) = varargin;
% Place optional args in memorable variable names
[nnsim, isAPoachEfficient, isDisplayOn, maxKAba, maxrAba, KSO, rSO] = optargs{:};
% ------- end checking input arguments.

disp('-|| Simulation procedure:')
disp('-|| Can sea otters and northern abalone co-exist under best management?')
disp('-|| Chades et al (under review)')

global DIR_results
global gfilename IS_DISPLAYED_GRAPH
global PARAM_LINEAR_FR

rand('state', 300);   % Same random seed

IS_DISPLAYED_GRAPH=isDisplayOn; % display 1= yes; 0=No

% Initialising/loading parameters
disp('-> Loading parameters');
load_param(isAPoachEfficient,maxKAba,maxrAba,KSO,rSO);
disp('<- Parameters loaded');

%% #File management
output_filename=[DIR_results,date,'_results_']; % For example: 02-Jun-2011_results_

% Local variables
nExpert=4;  % 4 experts
nFR=3;      % 3 functional responses: Sigmoid, Hyperbolic, Linear

for iExpert=1:nExpert % For all experts
    expert_name=['Expert',num2str(iExpert)];
    policy_iExpert = compute_expert_policy(iExpert);
    gfilename=fopen([output_filename,expert_name,'.txt'],'w+');
    fprintf(gfilename,'%d %f %f %f %f\n',isAPoachEfficient,maxKAba,...
        maxrAba,KSO,rSO);  % saving some common parameters with the file
    fprintf(gfilename,'%s \n| Average rew | +- | SO abundance | +- | Aba density | +- | Aba100mm | +- ', expert_name);
    disp(['<START>EXPERT ',num2str(iExpert)]);
    for iFR=1:nFR  % For all function response (FR)
        if iFR==1       % Sigmoid FR
            beginTag='Sig_';
            for iSig=1:2     % only 2 sig
                debutt=[beginTag,num2str(iSig)];
                disp(['<START> Simulation procedure for ',debutt]);
                 if IS_DISPLAYED_GRAPH==1, evaluate_strategy(policy_iExpert,[debutt,expert_name],iFR,iSig); end
                evaluate_strategy_mean(policy_iExpert,[debutt,expert_name],iFR,iSig,nnsim);
                disp(['<\END> Simulation procedure for ',debutt]);
            end
        elseif iFR==2   % Hyperbolic FR
            beginTag='Hyp_';
            for iHyp=1:2
                debutt=[beginTag,num2str(iHyp)];
                disp(['<START> Simulation procedure for ',debutt]);
                 if IS_DISPLAYED_GRAPH==1, evaluate_strategy(policy_iExpert,[debutt,expert_name],iFR,iHyp); end
                evaluate_strategy_mean(policy_iExpert,[debutt,expert_name],iFR,iHyp,nnsim);
                disp(['<\END> Simulation procedure for ',debutt]);
            end;
        elseif iFR==3   % Linear FR
            beginTag='Linear_';
            for iLinear=1:length(PARAM_LINEAR_FR)
                debutt=[beginTag,num2str(iLinear)];
                disp(['<START> Simulation procedure for ',debutt]);
                 if IS_DISPLAYED_GRAPH==1, evaluate_strategy(policy_iExpert,[debutt,expert_name],iFR,iLinear); end
                evaluate_strategy_mean(policy_iExpert,[debutt,expert_name],iFR,iLinear,nnsim);
                disp(['<\END> Simulation procedure for ',debutt]);
            end
        end
    end
    fclose(gfilename);
    disp(['Performance results saved in ',output_filename,expert_name,'.txt'])
    disp(['<END>EXPERT ',num2str(iExpert)]);
    close all
end
disp('-|| END Simulation procedure')
end
