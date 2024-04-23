% function
% [AbaPop,TAbaAdults,TSoPop,next_state,outcome,oilSpill]=simulation_t(AbaPop,SoPop,action,iFR,iFRnum)
% Action:
%   Simulation of conservation action (a) for t years on both populations
%   (AbaPop,SoPop) given a functional response (iFR,iFRnum).
% Input:
%   AbaPop: population vector of Abalone
%   SoPop:  sea otter abundance
%   action: management action
%   iFR,iFRnum: functional response parameters
% Output:
%   AbaPop: population vector of abalone
%   TAbaAdults: population vector of abalone over t time step.
%   TSoPop: sea otter abundance over t time step
%   next_state: corresponding MDP state
%   outcome: instant reward over t time step  
%   oilSpill:   information on oil spill over t time step
% Side effect:
%   This function was originally designed to accomodate a decision for t
%   time steps - useful for reinforcement procedure. Here I
%   assume t=1.
% Author: iadine.chades@csiro.au

function [AbaPop,AbaPopF,TAbaAdults,TSoPop,next_state,outcome,oilSpill]=simulation_t(AbaPop,AbaPopF,SoPop,action,iFR,iFRnum)
global PARAM_ABALONE  PARAM_QL 

t=PARAM_QL.t;
oilSpill=zeros(t,1);
TSoPop=zeros(t,1);
TAbaAdults=zeros(t,1);
outcome=zeros(t,1);

for i=1:t
    if SoPop<0
        error('We owe a big deal of otters')
    end;
    [UnUsed,SoPop,oil,AbaPop,AbaPopF]=perform_action(action,AbaPop,AbaPopF,SoPop,iFR,iFRnum);
    if SoPop<0
        error('We owe a big deal of otters') 
    end;
    oilSpill(i)=oil;
    TSoPop(i)=SoPop;
    TAbaAdults(i)=sum(AbaPop);
    AbaAdults=sum(AbaPop);
    next_state=[abaloneDensity2state(AbaAdults/PARAM_ABALONE.area), SOabundance2state(SoPop)];
    outcome(i)=Xreward(SoPop,AbaAdults/PARAM_ABALONE.area);
end


function reward=Xreward(so_pop,aba_pop)
% Joint rewards

global JOINT R PARAM_SO

K=PARAM_SO.k;

so_i=so_pop/K;
r=find(R(:,2)<=so_i);
r_so=R(r(end),3);

r=find(R(:,1)<=aba_pop);
r_aba=R(r(end),3);
if JOINT==1
    if r_so>r_aba
        reward=r_aba;
    else % r_so <= r_aba
        reward= r_so;
    end;
else
    reward=r_aba+r_so;
end