% function [cmin,cmax,dmin,dmax]=derive_sig_FR(k_aba_max)
% Action:
%    Fits 3 sigmoid FRs to observations/data.
%
% Input: 
%   k_aba_max: carrying capacity for abalone, 
%   Local variables:
%       a max number of prey eaten per day per sea otter (Pemax),
%       diving efficiency (effmin, effmax), 
%       a handling time (Th), and predation time (Tp).
% Output:
%   cmin,cmax,dmin,dmax: for best representation of uncertainty surrounding
%   the FRs we provide the min and max sigmoid FR (). These parameters
%   are then used to define and simulate the interaction between both
%   species.
% Side effects:
%   These functions are highly dependent on the interactions between 
%   species. We recommend defining your own functions if the application
%   varies.
%
% Author: iadine.chades@csiro.au
%

function [cmin,cmax,dmin,dmax]=derive_sig_FR(k_aba_max)
global IS_DISPLAYED_GRAPH

Pemax   = 18;
effmin  = 0.13;     
effmax  = 0.36;
Nmax=k_aba_max;
Tp=[0.25*effmin, 0.25*effmax, 0.25];
Th=111.75/3600/24;   
if IS_DISPLAYED_GRAPH==1
    figure('color','white','name','Sigmoid functional responses');
end
coul=[0.7 0.1 0];
N=0:0.01:Nmax;
Z=zeros(length(Tp),length(N));
Y=zeros(size(N));
c_sig=[];
d_sig=[];
theta=2;
for k=1:length(Tp)
    c=Tp(k)/Th;
    d2=c*Nmax^theta/Pemax-Nmax^theta;
    c_sig= [c_sig c];
    d_sig= [d_sig d2];
    for i=1:length(N)
        Z(k,i)= c*N(i)^theta/(d2+N(i)^theta);
        Y(i)=c*N(i)^theta/(d2+N(i)^theta);
    end
    if IS_DISPLAYED_GRAPH==1 plot(N,Y,'color',coul); 
    hold on;end
    coul=coul+[0.00 0.4 0.2];
end
if IS_DISPLAYED_GRAPH==1
    xlabel('Abalone density (m^-²)');
    ylabel('#Abalone eaten/day/sea otter');
    box off
    hold on
    plot(N,N*Pemax/Nmax,'--','color','k');
end
cmin=c_sig(1);
cmax=c_sig(3);
dmin=d_sig(1);
dmax=d_sig(3);
