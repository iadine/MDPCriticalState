% function [cost,nbso,oil,AbaPop,AbaPopF]=
% perform_action(a,AbaPop,AbaPopF,SoPop,iFR,iFRnum)
% Action:
%   Simulates the management decisions, NA growth, SO growth, Predation,
%   Poaching. Following the order of events :
%       0) Initialize conservation decision based on action (a)
%       1) Abalone population growth
%       2) SO population growth, culling and oil spill event
%       3) Abalone predation based on FR (iFR, iFRnum)
%       4) Abalone poaching
% Input:
%   a: management action,
%   AbaPop: abalone population vector,
%   AbaPopF: female abalone population vector,
%   SoPop: sea otter abundance
%   iFR,iFRnum: FR parameters
% Output:
%   cost: unused at this stage but could be included in future versions
%   nbso: updated so population
%   oil: oilspill
%   AbaPop: population vector updated.
% Side effect:
%   This procedure is the core of the simulation process. If errors occur
%   check the functions called.
% Author:   iadine.chades@csiro.au

function [cost,nbso,oil,AbaPop,AbaPopF]= perform_action(a,AbaPop,AbaPopF,SoPop,iFR,iFRnum)
global PARAM_SO
global PARAM_ABALONE

oil=0;
cost=0;     % No cost in this version.
nbso=SoPop;
poach=1;

if nbso <0
    error('Incorrect nbso value <0')
end

%% 0 - conservation decision
if a==0 || a==-1    % do nothing
    poach=1;
else
    if a==1 % introduction of Sea Otter only => impact_poaching = weak
        nbso=PARAM_SO.init_pop;
        poach=0; % WHY IS THIS 0 and not 1? 0 = optimistic case
    elseif a==2 % Anti-poaching only => impact_poaching = weak
        poach=0;
    elseif a==3 % Culling only => impact_poaching = normal
        poach=1;
    elseif a==4 % Half culling  - half antipoaching
        poach=0.5;
    end
end

%% 1 - growth of Aba pop/ density dependence
[AbaPop, AbaPopF] = northern_abalone_growth_t(AbaPopF);

%% 2 - SO
if nbso ~=0
    [nbso,oil]=sea_otter_growth(nbso);
    if a==3
        nbso= min(PARAM_SO.kculling*PARAM_SO.k,nbso); % (80% of K)
    end;
    if a==4 % 1/2s culling
        remove=(nbso-PARAM_SO.kculling*PARAM_SO.k)/2;
        nbso=min(nbso-remove,nbso);
    end
end
nbso=round(nbso);

%% 3 - predation
if nbso ~= 0
    AbaPop=predation_FR(nbso,AbaPop,iFR,iFRnum);
    AbaPopF= AbaPop*PARAM_ABALONE.rmf;
    if ~(AbaPop>0)
        AbaPop
        error('We are in debt, predators are starving');
    end
end

%% 4 - poaching
[AbaPop,AbaPopF]=compute_poaching_impact(AbaPop,AbaPopF,poach);
%AbaPop=round(AbaPop);

if ~(AbaPop>0)
    AbaPop
    error('We are in debt, we owe some northern abalone to poachers');
end



