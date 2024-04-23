% function [Ni,Nif]=northern_abalone_growth_t(Nif)
% Action:
%   Northern Abalone population model and dynamic
% Input:
%   Nif: female population vector (1..10)
% Output:
%   Ni: updated population vector
%   Nif: updated female population vector
% Global variables:
%   G: varies over time.
% Side effect:
%   reinitialise the matrix population before each simulation using
%   function initialising_northern_abalone
% Author:
%   iadine.chades@csiro.au

function [Ni1,Ni1f]=northern_abalone_growth_t(Ni1f) % todo
global PARAM_ABALONE G Ginit
global rand_k

% Initialising
area=PARAM_ABALONE.area;
k_aba_f=PARAM_ABALONE.kf;
r=PARAM_ABALONE.r;
k=k_aba_f*area;

rmax= max(eig(Ginit));
rmax=rand_k*r';
% Simulation
% Comment from here if you do not want stochastic growth rate.
x=rand;
if x<=rand_k(1)  % One way of simulating stochasticity on *r*
    rmax=r(1);
elseif x<=rand_k(1)+rand_k(2)
    rmax=r(2);
elseif x<=rand_k(3)+rand_k(1)+rand_k(2)
    rmax=r(3);
else
    rmax=r(4);
end
% end comment

% Density dependence
% all1=sum(Ni1f(6:10));
% r1=rmax*k*0.27/(rmax*all1-all1+k*0.27);
% max_eig_G=max(eig(G));
% m1=r1/max_eig_G;
% G(6:10,:)=G(6:10,:)*m1;
% essai Iadine/ JB

% k=k;
% %alld=sum(Ni1f(6:10));
all1=sum(Ni1f);
r1=rmax*k/(rmax*all1-all1+k);
max_eig_G=max(eig(G));
m1=r1/max_eig_G;
%G(6:10,:)=G(6:10,:)*m1;
G=G*m1;
% Growth
Ni1f=round(G*Ni1f);  % New population

% if sum(Ni1f(6:10))>(0.27*sum(Ni1f))
%     drate=(sum(Ni1f(6:10))/sum(Ni1f) - 0.27);
%     Ni1f(6:10)=(1-drate)*Ni1f(6:10);
%    % Ni1f(1:5)=(drate+1)*Ni1f(1:5);
% end
Ni1=Ni1f/PARAM_ABALONE.rmf;
end