% function [cmin,cmax,dmin,dmax]=derive_hyp_FR(k_aba_max)
% Action:
%    Fits 3 hyperbolic FRs to observations/data.
%
% Input: 
%   k_aba_max: carrying capacity for abalone, 
%   Local variables:
%       a max number of prey eaten (Pemax),
%       diving efficiency (effmin, effmax), 
%       a handling time (Th), and predation time (Tp).
% Output:
%   cmin,cmax,dmin,dmax: for best representation of uncertainty surrounding
%   the FRs we provide the min and max hyperbolic FR (). These parameters
%   are then used to define and simulate the interaction between both
%   species.
% Side effects:
%   These functions are highly dependent on the interactions between 
%   species. We recommend defining your own functions if the application
%   varies.
%
% Author: iadine.chades@csiro.au
%
function [cmin,cmax,dmin,dmax]=derive_hyp_FR(k_aba_max)
global IS_DISPLAYED_GRAPH

effmin=0.13;
effmax=0.36;
Pemax=18;       
Tp=[0.4*effmin, 0.4*effmax, 0.4*1];
Th=166/3600/24;

Nmax=k_aba_max;
if IS_DISPLAYED_GRAPH==1
    figure('color','white','name','Hyperbolic functional responses');
end
coul=[0.5 0.1 0];
N=0:0.01:Nmax;
Z=zeros(length(Tp),length(N));
Y=zeros(size(N));
c_hyp=[];
d_hyp=[];

for k=1:length(Tp)
    j=k;
    a=Pemax/(Nmax*Tp(k)-Pemax*Th*Nmax);
    
    c=Tp(k)/Th;
    d=1/(a*Th);
    c_hyp=[c_hyp c];
    d_hyp=[d_hyp d];
    for i=1:length(N)
        Z(k,i)= a*N(i)*Tp(k)/(1+a*Th*N(i));
        Y(i)=a*N(i)*Tp(k)/(1+a*Th*N(i));
    end
    if IS_DISPLAYED_GRAPH==1  plot(N,Y,'color',coul); 
    hold on;end
    coul=coul+[-0.2 0.3 0.5];
    %   end
end
if IS_DISPLAYED_GRAPH==1
    xlabel('Abalone density(m^-²)');
    ylabel('#Abalone eaten/day/sea otter');
    box off
    hold on
    plot(N,N*Pemax/Nmax,'--','color','k');
end
cmin=c_hyp(1);
dmin=d_hyp(1);
cmax=c_hyp(2);
dmax=d_hyp(2);
