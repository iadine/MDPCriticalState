% function function show_FR(iFR,iFRnum)
% Action: 
%   Plot the desired functional response (iFR,iFRnum)
% Input: 
%   iFR, iFRnum: Functional Response parameters
% Side effect:
%   Need to call figure() before calling show_FR
% Author: iadine.chades@csiro.au

function show_FR(iFR,iFRnum)
global PARAM_ABALONE PARAM_HYP_FR PARAM_SIG_FR PARAM_LINEAR_FR

k_max=PARAM_ABALONE.k;

N=0:0.01:k_max;
Z=zeros(1,size(N,2));

if iFR == 2 % Hyperbolic
    c=PARAM_HYP_FR.c(iFRnum);
    d=PARAM_HYP_FR.d(iFRnum);
    for i=1:size(N,2)
        Z(i)=c*N(i)/(d+N(i));
    end;
    plot(N,Z,'-','color','k','Linewidth',1);
    hold on;
%   text(0.8,4,'Hyperbolic FR f(N)=cN/(d+N)');
elseif iFR==1   % Sigmoid
    c=PARAM_SIG_FR.c(iFRnum);
    d=PARAM_SIG_FR.d(iFRnum);
    for i=1:size(N,2)
        Z(i)=c*N(i)^2/(d+N(i)^2);
    end;
    plot(N,Z,'-','color','k','Linewidth',1);
    hold on;
    %text(0.8,4,'Sigmoid FR f(N)=cN²/(d+N²)');
elseif iFR==3       % Linear
    d_max=PARAM_LINEAR_FR(iFRnum,1);
    Pemax=PARAM_LINEAR_FR(iFRnum,2);
    for i=1:size(N,2)
        if N(i)<d_max
        Z(i)=Pemax*N(i)/d_max;
        else
            Z(i)=Pemax;
        end
    end;
    plot(N,Z,'-','color','k','Linewidth',1);
    hold on;
end;
xlabel('Abalone density (m^{-2})')
ylabel('No. abalone eaten/Sea otter/Day');
hold on
%plot(N,N*18/k_max,'-','color','k');
xlim([0 k_max])
box off
