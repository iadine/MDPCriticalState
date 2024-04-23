% function Tr=ComputeTransition(iFR,iFRnum,ksim)
% Action:
%   Computes the transition probability matrices by simulation that will be
%   used to compute the optimal or near optimal strategy.
% Input:
%   iFR: functional response family [1,3]
%   iFRnum: functional response number within its family
%   ksim:   number of simulations for each state action pair.
% Output:
%   Tr: matrix (nbs,nbs,PARAM_MDP.nb_action)
% Side effect:
%   The sea otter population might go extinct due to oil spills and message
%   might appear in the console. 
%   Because it is a simulation procedure this function can take a long
%   time. It is recommended to use small values of ksim for approximate but
%   fast results (ksim<100); and large number of ksim for more accurate
%   estimation (ksim>1000).
%
% Author: iadine.chades@csiro.au
%

function Tr=compute_transition(iFR,iFRnum,ksim)

global PARAM_MDP G PARAM_SO

Ginit=G;
nbsa=PARAM_MDP.nbs_aba;
nbso=PARAM_MDP.nbs_so;
nbs=PARAM_MDP.nbs_aba*PARAM_MDP.nbs_so;   % nbs = number of states
limit_state_culling=SOabundance2state(PARAM_SO.kculling*PARAM_SO.k); %first state we can do culling

%% initialize Tr-matrix
Tr=zeros(nbs,nbs,PARAM_MDP.nb_action);
%Rew=zeros(nbs,1);
for i=1:nbs
    s=seeState(i,PARAM_MDP.state_matrix);
    if s(2)~=0         % no otter
        Tr(i,:,1)=NaN; % forbid action 1 for state with SO pop
    end
    if s(2)<limit_state_culling
        Tr(i,:,3)=NaN; % forbid action 3 and 4 for state with SO pop below limit_state_culling
        Tr(i,:,4)=NaN;
    end;
end;

for i=0:nbsa-1
    % i current abalone state
    for j=0:nbso-1
        %j current so state
        if j==0
            TAction=[1,2];
        elseif j<limit_state_culling
            TAction=[2];
        elseif j>=limit_state_culling
            TAction=[2,3,4];
        end
        % for each action, simulate the transition
        for a=1:length(TAction)
                for k=1:ksim
                    % random pick of an abalone value in this state
                    G=Ginit;
                    [aba_density,AbaPop, AbaPopF]= randnumberaba(i);
                    % random pick of a sea otter abundance in this state
                    nb_so= randnumberso(j);
                    Init_state= seeIndex([abaloneDensity2state(aba_density),SOabundance2state(nb_so)]);                   
                    [unUsed0,unUsed1,unUsed2,unUsed3,next_state,unUsed4,unUsed5]=simulation_t(AbaPop,AbaPopF,nb_so,TAction(a),iFR,iFRnum);
                    Tr(Init_state,seeIndex(next_state),TAction(a))=Tr(Init_state,seeIndex(next_state),TAction(a))+1/ksim;
                end
        end
    end
end
G=Ginit;
%save 'temp_tr.mat' Tr
end

function [aba_density,AbaPop,AbaPopF]=randnumberaba(statea)
global PARAM_ABALONE

min_dens=state2dens_aba(statea);
max_dens=state2dens_aba(statea+1);

aba_density=(max_dens-min_dens)*rand+min_dens;
if aba_density<=0.001
    aba_density=0.001;
elseif aba_density>=PARAM_ABALONE.k
    aba_density
end
if abaloneDensity2state(aba_density)~=statea
    'mismatch! bad conversion'
end
[AbaPop, AbaPopF]=initialising_northern_abalone(0,aba_density); %update G sinon ca part en cacahuete dans la simu!
end

function nb_so=randnumberso(stateso)
min_so=state2abund_so(stateso);
max_so=state2abund_so(stateso+1);
if min_so~=-1
    nb_so=(max_so-min_so)*rand+min_so;
else
    nb_so=0;
end
if SOabundance2state(nb_so)~=stateso
    'mismatch! mauvaise conversion so'
end
end

function densitymin=state2dens_aba(s)
global PARAM_ABALONE PARAM_MDP

densitymin=s*PARAM_ABALONE.discUnit;
if s==PARAM_MDP.nbs_aba
   densitymin=PARAM_ABALONE.k;
end
end

function min_so=state2abund_so(s)
global PARAM_MDP PARAM_SO

K=PARAM_SO.k;
nbs=PARAM_MDP.nbs_so;
tranche=K/(nbs-1);
if s==0
    min_so=-1;
else
    min_so=(s-1)*tranche;
end
end