% Function new_Tabundance_prey=
% predation_FR(abundance_predator,Tabundance_prey,iFR, iFRnum)
% Action:
%   Computes the new abundance of the prey population after predation 
%   following a functional response type (iFR) and number (iFRnum).
% Input:
%   abundance_predator: sea otter abundance
%   Tabundance_prey:    northern abalone population vector
%   iFR,iFRnum:     functional response
% Output:
%   new_Tabundance_prey: updated northern abalone population vector
% Side effect:
%   We assume a uniform predation across the population, however when a
%   life stage is depleted we remove preys from a different life-stage.
%   
% Author: iadine.chades@csiro.au

function new_Tabundance_prey = predation_FR(abundance_predator,Tabundance_prey,iFR,iFRnum)
global PARAM_ABALONE PARAM_LINEAR_FR PARAM_SIG_FR PARAM_HYP_FR

area_aba=PARAM_ABALONE.area;
if abundance_predator == 0    % no predator
    new_Tabundance_prey=Tabundance_prey;
else            % predator
    days=365;   % 1 year = 365 days
    sum_Tabundance_prey=sum(Tabundance_prey);
    Nd=sum_Tabundance_prey/area_aba;  % Nd = density
   if iFR==2 %hyp
        c=PARAM_HYP_FR.c(iFRnum);
        d=PARAM_HYP_FR.d(iFRnum);
        removed_prey=c*Nd/(d+Nd)*days*abundance_predator;
    elseif iFR==1 % sigmoid
        c=PARAM_SIG_FR.c(iFRnum);
        d=PARAM_SIG_FR.d(iFRnum);
        removed_prey=c*Nd^2/(d+Nd^2)*days*abundance_predator; % good
    elseif (iFR==3) || (iFR==0) %linear 
        d_max=PARAM_LINEAR_FR(iFRnum,1);
        Pemax=PARAM_LINEAR_FR(iFRnum,2);
        if Nd<d_max 
            removed_prey=Pemax*Nd/d_max*days*abundance_predator;
        else removed_prey=Pemax*Nd/d_max*days*abundance_predator;
        end
    else    % error
        disp('Incorrect FR type');
        error('FR type incorrectly defined');
    end;
    if removed_prey>sum_Tabundance_prey
        'population crash predationFR(sig,abundance_predator,Tabundance_prey,area_aba,c,d)'
        new_Tabundance_prey=zeros(size(Tabundance_prey,1),1);
    else
        removed=ones(size(Tabundance_prey,1),1)*removed_prey/10; % assume uniform predation across prey population
        new_Tabundance_prey=Tabundance_prey-removed;
        if (sum(new_Tabundance_prey<0)>0)   % check if population is <0
            i=1;
            while (sum(new_Tabundance_prey<0)>0)
                if new_Tabundance_prey(i)<0
                    index=mod(i+1,size(new_Tabundance_prey,1))+1;
                    new_Tabundance_prey(index)= ...
                        new_Tabundance_prey(index)-new_Tabundance_prey(i);
                    new_Tabundance_prey(i)=0;
                end;
                i=i+1;
                if i>size(new_Tabundance_prey,1)
                    i=1;
                end;
            end
        end
    end
end