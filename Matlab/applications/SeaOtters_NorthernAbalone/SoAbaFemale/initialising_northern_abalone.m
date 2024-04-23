% function N1=initialising_northern_abalone(disp,target_aba_density)
% Action:
%   Initialisation of the northern abalone matrix population model. This function
%   reset the population model to initial values and stabilise the
%   population before simulations are performed. Because our model is
%   density dependent G changes over time and must be reinitialised.
% Input:
%   'disp' calls the plot function
%   target_aba_density: used as stable density of abalone to reach after T 
%   init time step (male and female).
% 
% Output:
%   N1: a population vector of (1..10), represents both male and female
%   N1f: a FEMALE population vector of (1..10).
% Global:
%   growth rate r=PARAM_ABALONE.r;
%   init time step T=PARAM_ABALONE.T_init;
%   
% Side effect:
%   The important parameters are defined in load_parameter(). Please refer
%   to this function. If Ni1 is changed and updated for a different
%   application please do change the local values defined here.

function [N1,N1f]=initialising_northern_abalone(disp,target_aba_density)
global PARAM_ABALONE Ginit G
global rand_k

area=PARAM_ABALONE.area;
ratio_mf=PARAM_ABALONE.rmf;
if exist('target_aba_density','var')
    k_target=target_aba_density*ratio_mf;
else
    k_target=PARAM_ABALONE.init_pop*ratio_mf;
end

r=PARAM_ABALONE.r;
T=PARAM_ABALONE.T_init;

% Because G is dynamic we need to reinitialise its value
% Transition matrix initialization
G= Ginit;

% Initial population - could be defined as global
N1f=[0.047;0.056;0.040;0.023;0.018;0.007;0.011;0.003;0.00;0.025]*area*ratio_mf;    %Female only

% Abapop stores the population over time
AbaPop=zeros(4,T);
k=k_target*area;
rmax=rand_k*r';

% Simulation over T years to stabilise population
for i=1:T
    AbaPop(1,i)=sum(N1f(1:10));  % all
    all1=AbaPop(1,i);
    r1=rmax*k/(rmax*all1-all1+k);
    m1=r1/max(eig(G));
    G=G*m1;
    N1f=G*N1f;
end

N1=N1f/ratio_mf;    % N1 represents both male and female population
% Plot initialisation if asked for
if disp==1
    figure
    plot((AbaPop(1,:)/area)','c-x');
    xlabel('time (yrs)')
    ylabel('density (m^2)');
    legend('scramble (R)','Location', 'best');
end
