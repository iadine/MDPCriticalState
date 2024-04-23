% function main_SDP_SA(nsim,                % number of simulation
%                   isAPoachEfficient,   % 0 or 1    reduce poaching by 0.5/0.75 (default 0)
%                   isDisplayOn,         % 0 or 1    figures off/on   (default 0)
%                   maxKAba,             % carrying capacity in best habitat (default 3.34)
%                   maxrAba,             % max growth rate in best habitat (default 1.6)
%                   KSO,                 % carrying capacity sea otter (default 4073)
%                   rSO)                 % growth rate sea otter (default 0.191)
%
% main_SDP_SA is the function to run the sensitivity analysis on the optimisation procedure.
% main_SDP_SA(1000,0,0,3.34,1.6,4073,0.191)
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
% - The simulation results are stored in Results/[date]_results_SA_VI.txt
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

function main_SDP_SA(varargin)
rand('state', 300);     % Same random seed

%------- checking input arguments
numvarargs = length(varargin);
if numvarargs > 7
    error('myfuns:main_SDP:TooManyInputs', ...
        'requires at most 7 optional inputs');
end
% set defaults for optional inputs
optargs = {500 0 0 3.34 1.6 4073 0.191};
optargs(1:numvarargs) = varargin;
% Place optional args in memorable variable names
[nsim, isAPoachEfficient, isDisplayOn, maxKAba, maxrAba, KSO, rSO] = optargs{:};
% ------- end checking input arguments.

disp('-|| Sensitivity analusis of the optimisation procedure:')
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
output_filename=[DIR_results,date,'_results_SA_VI']; % For example: 02-Jun-2011_results_SA_VI
gfilename=fopen([output_filename,'.txt'],'w+');
fprintf(gfilename,'%d %f %f %f %f\n',isAPoachEfficient,maxKAba,...
    maxrAba,KSO,rSO);  % saving some common parameters with the file
fprintf(gfilename,'SA_VI \n| Average rew | +- | SO abundance | +- | Aba density | +- | Aba100mm | +- ');
% End file management

nFR=3;  % There are 3 families of functional responses (FR):
% sigmoid, hyperbolic and linear.

new_k=[3.34-0.2,3.34 , 3.34+0.2];
new_r=[1.6-0.05, 1.6, 1.6+0.05];
for ik=1:length(new_k)
    for ir=1:length(new_r)
        disp('-> Loading parameters');
        load_param(isAPoachEfficient,new_k(ik),new_r(ir),KSO,rSO);
        disp('<- Parameters loaded');
        for iFR=3:nFR
            if iFR==3   % Linear FR ONLY for sensitivity analysis!
                beginTag='Linear_';
                for iLinear=1:length(PARAM_LINEAR_FR) % variable number of Linear FR
                    debutt=[beginTag,num2str(iLinear)];
                    disp(['<START> Optimisation procedure for functional response ',debutt]);
                    Tr= compute_transition(iFR,iLinear,ksim);   % Estimation of Tr
                    Policy_VI= value_iteration();               % Optimisation algorithm
                     evaluate_strategy_mean(Policy_VI,[debutt,'VI'],iFR,iLinear,k_max);
                    disp(['<\END> Optimisation procedure for ',debutt]);
                end
            end
        end
    end
end
fclose(gfilename);
disp(['Performance results saved in ', output_filename,'.txt'])
disp('-|| END Optimisation procedure')

end