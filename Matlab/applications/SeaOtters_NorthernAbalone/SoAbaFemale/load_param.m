% function load_param( isAPoachEfficient,   % 0 or 1    reduce poaching by 0.5/0.75 (default 0)
%                      maxKAba,             % carrying capacity in best habitat (default 3.34)
%                      maxrAba,             % max growth rate in best habitat (default 1.6)
%                      KSO,                 % carrying capacity sea otter (default 4073)
%                      rSO)                 % growth rate sea otter (default 0.191)
% Action:
%   initialize most ot the important parameters for the population models, 
%   functional response, optimisation models, joint reward function.
% Input:
%   low = 1; poach_low=0.1; 10% poaching intensity when anti-poaching implemented 
%   low ~= 1; poach_low=0.05; % 5% poaching intensity when anti-poaching 
% Output:
%   Global variables: PARAM_ABALONE PARAM_LINEAR_FR PARAM_SIG_FR 
%   PARAM_HYP_FR PARAM_POACHING PARAM_SO PARAM_MDP PARAM_QL R TOTAL 
%   DIR_results DIR_fig Ginit Ginitl Ginith G Gl Gh rand_k JOINT
%   All these parameters are defined according to the litterature or data
%   collected unless stated otherwise. Please refer to the paper Chades et
%   al (in review) for details.
% Side effect:  
%   This program is specifically designed for the sea otter and abalone
%   problem. Consider contacting the author if you would like to use this
%   program for a different problem.
%   Most global variables can be changed to assess the sensitivity of the
%   results. In particular: i) the abalone and sea otter population
%   parameters, ii) the functional response iii) the objective/reward 
%   function of the optimisation procedure.
%   However our results suggest that poaching activity is the main threat
%   to abalone and removing sea otters has a marginal impact on the density
%   of abalone. 
%
% Author: iadine.chades@csiro.au

function load_param(isAPoachEfficient,maxKAba,maxrAba,KSO,rSO)                

global PARAM_ABALONE PARAM_LINEAR_FR PARAM_SIG_FR PARAM_HYP_FR PARAM_POACHING PARAM_SO PARAM_MDP 
global PARAM_QL R TOTAL DIR_results DIR_fig INIT_N

DIR_fig='Fig/'; if ~exist(DIR_fig,'dir') mkdir(DIR_fig); end
DIR_results='Results/'; if ~exist(DIR_results,'dir') mkdir(DIR_results); end


% ratio male and female abalone
ratio_mf=0.5;

% parameters that can be change
TOTAL=100000;       % the number of learning trials for the Reinforcement Learning
T=100;              %  length of a managing episode (usually 100 years)

% parameters that are not recommended to change
t=1;                % t is the time step at which a decision is made  REM= obsolete

%% #Abalone population variables
% Growth rate of 4 populations from bad quality habitat to highly suitable
rmax=derive_rmax(maxrAba);      % rmax=[1.05,1.2,1.4,1.6];  
k_aba_h= (derive_kabah(maxKAba)); % k_aba_h=[0.837,1.67,2.5,3.34];  

% do not change/ specific to pacific rim national park, BC, Canada.
surf_k=[64.78 109.47 90.10 31.35]*(10^6); % surfaces of the 4 kinds of habitat
area_aba=sum(surf_k);                   % total area of abalone

k_aba=sum(k_aba_h.*surf_k)/area_aba;    % Average carrying capacity
k_aba_fem= k_aba*ratio_mf;         % Average carrying capacity
k_init=0.23;                       % To initialize population
T_init=50;                         % Time to initialize population

% management objective 
global JOINT
JOINT = 1;
global rand_k
rand_k=surf_k./area_aba;
global Ginit Ginitl Ginith
global G Gl Gh

% Initial population 
INIT_N=[0.047;0.056;0.040;0.023;0.018;0.007;0.011;0.003;0.00;0.025]*area_aba*ratio_mf;    %Female only

% Female abalone population matrix initialization
finit=[0.136;0.26;0.38;0.491;0.593;0.683;0.745;0.795;0.835;1.166];
finit=finit*(10^6)*ratio_mf; % female only
sk=0.818;
sj=[5.42*10^-7, 5.42*10^-7+3.15*10^-7, 5.42*10^-7-3.15*10^-7];

sk=sk*(1+0);              % taking of estimated fishing pressure
sj=sj*(1+0);              % Taking of estimated fishing pressure 

sj_m=sj(1);                 % Average survival rate
sj_l=sj(3);                 % Low survival rate
sj_h=sj(2);                 % High survival rate




G=zeros(10);
G(1,:)= sj_m.*finit;
G(2,1)= sk; G(3,2)=sk; G(4,3)=sk;  G(5,4) =sk; G(6,5)  =sk;
G(7,6)= sk; G(8,7)=sk; G(9,8)=sk;  G(10,9)=sk; G(10,10)=sk;
Ginit=G;

% Note: Using Gl and/or Gh do not make a significant difference
% Gl is G with low survival rate
Gl=zeros(10);
Gl(1,:)= sj_l.*finit;
Gl(2,1)= sk; Gl(3,2)=sk; Gl(4,3)=sk;  Gl(5,4) =sk; Gl(6,5)  =sk;
Gl(7,6)= sk; Gl(8,7)=sk; Gl(9,8)=sk;  Gl(10,9)=sk; Gl(10,10)=sk;
Ginitl=Gl;
% Gh is G with high survival rate
Gh=zeros(10);
Gh(1,:)= sj_h.*finit;
Gh(2,1)= sk; Gh(3,2)=sk; Gh(4,3)=sk;  Gh(5,4) =sk; Gh(6,5)  =sk;
Gh(7,6)= sk; Gh(8,7)=sk; Gh(9,8)=sk;  Gh(10,9)=sk; Gh(10,10)=sk;
Ginith=Gh;

Ginit=Gh;

% Abalone population system states for optimisation model
discretisationUnit=0.05; % discretisation step at which the decision is made
nbStates_abalone=ceil(k_aba/discretisationUnit);  
         % number of states representing the status of the species
         % # of states representing the Abalone population nbs_aba=|Xa|                                
PARAM_ABALONE=struct('k',k_aba,'kf',k_aba_fem,'r',rmax,'init_pop',k_init,'T_init',T_init,...
    'area',area_aba,'poaching',1,'discUnit',discretisationUnit,'rmf',ratio_mf);

%% #poaching parameters on Abalone
% can be changed
if isAPoachEfficient==0
    poach_low=0.1;  % decrease in 50% of poaching intensity when anti-poaching
else
    poach_low=0.05; % decrease in 75% of poaching intensity when anti-poaching
end
% do not change
% current estimate of poaching activity
poach_thr=0.01;  % Poaching threshold e.g. density of abalone is very low 
poach_high=0.23; % Poaching intensity 0.20+-0.3 
poach_med=0.17;  % Poaching intensity 0.20+-0.3
PARAM_POACHING=struct('thr',poach_thr,'high',poach_high,'med',poach_med,...
                      'low',poach_low);

%% #Sea otter population variables
k_so=KSO;               % 4073- Carrying capacity of the pacific rim national park (Gregr et al. 2008)
if k_so <500
    error('k_so too small program not designed for small populations');
end
area_so=1036*10^6;      % Estimated area of suitable habitat (Gregr et al. 2008)
r_so=rSO;               % 0.191- Growth rate based on (Nichol 2007)
init_so=100;            % Number of individuals initially introduced 
authorise_culling=0.6;  % Authorise culling when SO has reached 0.6*k_so
oil_spill_frequency=0.1;% Gerber et al (2004) 
dead_prct_min=0.23;     % Gerber et al (2004)
dead_prct_max=0.42;
extinct_when=10;    % threshold of extinct population (abundance)

PARAM_SO=struct('k',k_so,'r',r_so,'init_pop',init_so,'area',area_so,...
                'os_freq',oil_spill_frequency,'allee_effect',extinct_when,...
                'kculling',authorise_culling,'dead_prct_min',dead_prct_min,...
                'dead_prct_max',dead_prct_max);

%% #Functional response parameters
[c1_h,c2_h,d1_h,d2_h]= derive_hyp_FR(k_aba); % deriving the hyperbolic functional response based on k_aba
[c1_s,c2_s,d1_s,d2_s]= derive_sig_FR(k_aba); % deriving the sigmoid functional response based on k_aba
    
Pemax=18; % Maximum of prey eaten per day. Tinker et al (2006)
PARAM_SIG_FR=struct('c',[c1_s,c2_s],'d',[d1_s,d2_s],'Pemax',Pemax,...
                    'd_max',k_aba);
PARAM_HYP_FR=struct('c',[c1_h,c2_h],'d',[d1_h,d2_h],'Pemax',Pemax,...
                    'd_max',k_aba);

% set of 4 linear FRs defined based on Pemax and k_aba
PARAM_LINEAR_FR=[ 
    k_aba Pemax/3; 
    k_aba 2*Pemax/3; 
    2*k_aba/3 Pemax; 
    k_aba/3 Pemax;];


%% #Sea otter population system states
nbStates_so=21;  % # of states representing the SO population nbs_seaotter=|Xs|

%% #Conservation decisions
nba=4;  % # of actions {1=intro_so, 2=anti_poaching, 3=culling, 4= both} 

%% #Learning algorithms parameters
exploration_Q=0.9;  % Q-learning exploration rate 1-0.9
alpha_Q=0.1;        % Q-learning learning rate

klim=200;      % Simulation runs for ongoing evaluations
freq=1000;     % Frequence of ongoing evaluations

%% Build the state matrix from abundance and density to x_a,x_so
state_matrix=build_state_matrix(nbStates_abalone,nbStates_so);

PARAM_MDP=struct('nbs_aba',nbStates_abalone,'nbs_so',nbStates_so,...
                 'nb_action',nba,'state_matrix',state_matrix,'epsi',0.001);
PARAM_QL=struct('total',TOTAL,'exploration_Q',exploration_Q,...
                'alpha_Q',alpha_Q,'Time_Horizon',T,'klim',klim,...
                'freq',freq,'t',t);

%% reward function
% Do not change. The levels identify the important stages of the recovery
% process

level_aba=[0 0.02 0.1 0.5 k_aba]; % Female and Male
level_so= [0 0.3  0.4 0.6 1 ]; 
R=build_joint_reward(k_aba,level_so,level_aba);

display(PARAM_ABALONE) 
display(PARAM_SO) 
display(PARAM_POACHING) 
display(PARAM_LINEAR_FR)
display(PARAM_SIG_FR)
display(PARAM_HYP_FR)
display(PARAM_MDP)
display(PARAM_QL)
end

function rmax=derive_rmax(maxrAba)  % default value is rmax~[1.05,1.2,1.4,1.6];  
 if maxrAba<1.05
     error('maxrAba has to be > 1.05')
 end
b=1.05;             % growth rate in poor habitat
a= (maxrAba-b)/3;   % assumed a linear relationship between different habitat type
rmax=[b,a+b,2*a+b,maxrAba];
end

function k_aba_h= derive_kabah(maxKAba) % k_aba_h=[0.837,1.67,2.5,3.34];  % new K
a=maxKAba/4;  % assumed a linear relationship between different habitat type
k_aba_h=[a, 2*a, 3*a, maxKAba];
end

