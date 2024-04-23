% function N=compute_poaching_impact(N,poach) Compute the impact of illegal fishing activity (poach)
% Action:
%   Compute the impact of illegal fishing activity (poach) on the 
% Northern Abalone population (N).
% Input: 
%   N: abalone population vector (1:10).
%   Nf: female abalone population vector (1:10).
%   poach:  >0  => high and medium poaching
%           ==0 => low poaching.
% Output:
%   N: updated abalone population after poaching activity.
%   Nf: updated female abalone population after poaching activity.
% Side effect:
%   if the density of abalone has reached a threshold PARAM_POACHING.thr
%   defined in function param_load(); the illegal fishing activity is
%   assumed to have an impact of 0.01 (1%).
% Author: iadine chades iadine.chades@csiro.au
% Copyright: See README.txt

function [N,Nf]=compute_poaching_impact(N,Nf,poach)
global PARAM_POACHING PARAM_ABALONE

h=PARAM_POACHING.high;
m=PARAM_POACHING.med;
L=PARAM_POACHING.low;
area_aba=PARAM_ABALONE.area;
thr=PARAM_POACHING.thr;

all1=sum(N(1:10));
impact_poaching=rand;
if poach>0  % High and Medium poaching
    if impact_poaching>0.5
        impact=(h-L)*poach+L;
    else impact=(m-L)*poach+L;
    end;
elseif poach==0 % Low poaching
    if impact_poaching>0.5
        impact=L+0.01;
    else impact=L-0.01;
    end
end
if (all1/area_aba) < thr % Northern abalone density has reached threshold
    impact=0.01;
end
N=(1-impact).*N;
Nf=(1-impact).*Nf;
end