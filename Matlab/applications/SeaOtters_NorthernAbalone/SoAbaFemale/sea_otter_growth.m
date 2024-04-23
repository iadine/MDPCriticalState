% Function [Y,OS]=sea_otter_growth(N)
% Action:
%   Sea otter population model including risk of oil spill
%   Ricker model from Gerber et al (2004)
% INPUT: sea otter abundance (N)
% OUTPUT: sea otter abundance (Y), Intensity (OS)
%
% Author: iadine.chades@csiro.au
function [Y,OS]=sea_otter_growth(N)

global PARAM_SO
K=PARAM_SO.k;
r=PARAM_SO.r;
min_d=PARAM_SO.dead_prct_min;
max_d=PARAM_SO.dead_prct_max;
oil_spill_frequency=PARAM_SO.os_freq;

Theta=1; % change this parameter for theta logistic
OS=0;
Y=N*exp(r*(1-N/K)^Theta);
oil_spill=rand;
if oil_spill < oil_spill_frequency
    dead_prct=min_d+(max_d-min_d)*rand;
    OS=dead_prct;
    Y=Y-Y*dead_prct;
    if Y<PARAM_SO.allee_effect % allee effect
        disp('Population of sea otter goes extinct bloody oil spill :-[');
        Y=0;
    end
end;



