% function Pib=value_iteration()
% Action:
%   Performs the value iteration algorithm
% Input:
%   Requires a global variable Tr (computed from the function 
%   compute_transition); PARAM_SO PARAM_MDP
% Output:
%   Pib: the strategy epsilon-optimal state-> action
%   global variable: q
% Side effect: 
%   None. Consider using the MDP toolbox from INRA if you require a
%   different optimisation algorithm.
%   http://www.inra.fr/internet/Departements/MIA/T/MDPtoolbox/
%
% Author: iadine.chades@csiro.au

function Pib=value_iteration()
global PARAM_SO PARAM_MDP Tr q

epsi=PARAM_MDP.epsi;
nbs=PARAM_MDP.nbs_aba*PARAM_MDP.nbs_so;     % nbs = number of states
gamma=0.9;                                  % gamma = discount factor

limit_state_culling=SOabundance2state(PARAM_SO.kculling*PARAM_SO.k); %first state we can do culling

%% initialize Q-matrix
q=zeros(nbs,PARAM_MDP.nb_action);
Vt=zeros(nbs,1);
V=zeros(nbs,1);
for i=1:nbs
    s=seeState(i,PARAM_MDP.state_matrix);
    if s(2)~=0        % no otter
        q(i,1)=NaN; % forbid action 1 for state with SO pop
    end
    if s(2)<limit_state_culling
        q(i,3)=NaN; % forbid action 3 and 4 for state with SO pop below limit_state_culling
        q(i,4)=NaN;
    end;
end;

bidule=1;
while bidule
    for i=1:nbs
        if s(2)==0
            TAction=[1,2];
        elseif s(2)<limit_state_culling
            TAction=[2];
        elseif s(2)>=limit_state_culling
            TAction=[2,3,4];
        end
        for a=1:length(TAction)
            q(i,TAction(a))=XSReward(i)+gamma*sum(Tr(i,:,TAction(a))*Vt(:));
        end
        V(i) = max(q(i,:));
    end
    if (abs(V-Vt))<epsi
        bidule=0;
    end
    Vt=V;
end
Pib=zeros(nbs,1);
for j=1:nbs
    [h,Pib(j)]=max(q(j,:));
    if isnan(h)
        Pib(j)=-1;
    end
end
end

function reward=XSReward(index)
% Joint rewards
% initialize values
global JOINT R PARAM_SO PARAM_MDP
K=PARAM_SO.k;

s=seeState(index,PARAM_MDP.state_matrix);
statea=s(1); stateso=s(2);
min_dens=state2dens_aba(statea);
max_dens=state2dens_aba(statea+1);

aba_pop=min_dens+(max_dens-min_dens)/2;
if aba_pop<=0
    'We owe a big deal of abalone XSReward'    
    aba_pop=0.001;
end

min_so=state2abund_so(stateso);
max_so=state2abund_so(stateso+1);
if min_so~=-1
    so_pop=min_so+(max_so-min_so)/2;
else
    so_pop=0;
end
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