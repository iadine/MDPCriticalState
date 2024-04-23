%Quarantine, Surveillance, Control SDP

%States:
%State 1 = species is absent
%State 2 = species is localised and undetected
%State 3 = species is localised and detected
%State 4 = species is widespread (always detected)

%Parameters
T = 10; %Time horizon

CW = 2900000; %Impact of widespread population
k = 0.01; %Ratio of impact of widespread:localised population
CL = k*CW; %Impact of localised population
CE = 96667; %Cost of eradicating localised population
Pl = 1.0;

p0 = 0.99; %Probability of reinvasion if no quarantine
alpha = 0.00000207; %Effectiveness of quarantine
beta = 0.0000157; %Effectiveness of surveillance
lambda = 0.00000523; %Effectiveness of eradication when widespread
g = 0.5;

Trans = zeros(4, 4); %Transition matrix
TransCost = zeros(4, 4); %Transition matrix
OptDec = zeros(4, 3, T, 2); %Array to store optimal allocations to Q, S, C for each state in each time step
OptCost = zeros(4, T, 2); %Array to store net expected cost for each state in each time step
AbsDec = zeros(6, T); %Array to store results for absent state, in format that can be easily graphed

inc = 50000;
VecL = CW/inc;

%SDP
for j = 1:2 %two different discount rates, 0% and 10%
    
    r = (j-1)*0.1;
        
    FutureOptNEC = [0, CL/((1+r).^(T+1)), (CL+CE)/((1+r).^(T+1)), CW/((1+r).^(T+1))]; %Initialising future costs in time T+1 (discounted)
    
    %Loop over all time steps
    for t = T:-1:1
        
        disc = 1/((1+r).^t); %discount variable, multiply cost in this timestep by this to discount
        
        %Loop over all possible states at time t
        for n = 1:4
            OptNEC = inf;
            
            %Loop over all possible allocations to quarantine
            for q = 1:VecL+1
                
                %Loop over all possible allocations to surveillance
                for s = 1:VecL+1
                    
                    %Loop over all possible allocations to control
                    for c = 1:VecL+1
                        
                        SumECost = 0;
                        
                        %Allocations
                        Xq = (q-1)*inc;
                        Xs = (s-1)*inc;
                        Xc = (c-1)*inc;
                        
                        %Probabilities
                        Pi = p0*exp(-Xq*alpha);
                        Pd = 1 - exp(-Xs*beta);
                        Pw = 1 - exp(-Xc*lambda);
                        
                        for m = 1:4 %Loop over all possible states at time t+1
                            if n==1
                                
                                if m==1
                                    Trans(1, 1) = (1 - Pi) + (Pi*Pd*Pl) + (Pi*(1-Pd)*g*Pw);
                                    TransCost(1, 1) = (Pi*Pd*(CL+CE)*Pl) + (Pi*(1-Pd)*g*CW*Pw);
                                elseif m==2
                                    Trans(1, 2) = (Pi*(1-Pd)*(1-g)) + (Pi*Pd*(1-Pl));
                                    TransCost(1, 2) = (Pi*(1-Pd)*(1-g)*CL*disc) + (Pi*Pd*(CL+CE)*(1-Pl));
                                elseif m==3
                                    Trans(1, 3) = 0;
                                    TransCost(1, 3) = 0;
                                else %m==4
                                    Trans(1, 4) = Pi*(1-Pd)*g*(1-Pw);
                                    TransCost(1, 4) = Pi*(1-Pd)*g*CW*(1-Pw);
                                end
                                
                            elseif n==2
                                
                                if m==1
                                    Trans(2, 1) = (Pd*Pl) + ((1-Pd)*g*Pw);
                                    TransCost(2, 1) = (Pd*(CL+CE)*Pl) + ((1-Pd)*g*CW*Pw);
                                elseif m==2
                                    Trans(2, 2) = ((1-Pd)*(1-g)) + (Pd*(1-Pl));
                                    TransCost(2, 2) = ((1-Pd)*(1-g)*CL*disc) + (Pd*(CL+CE)*(1-Pl));
                                elseif m==3
                                    Trans(2, 3) = 0;
                                    TransCost(2, 3) = 0;
                                else %m==4
                                    Trans(2, 4) = (1 - Pd)*g*(1-Pw);
                                    TransCost(2, 4) = (1 - Pd)*g*CW*(1-Pw);
                                end
                                
                            elseif n==3
                                
                                if m==1
                                    Trans(3, 1) = Pl;
                                    TransCost(3, 1) = (CL+CE)*Pl;
                                elseif m==2
                                    Trans(3, 2) = (1-Pl);
                                    TransCost(3, 2) = (CL+CE)*(1-Pl);
                                elseif m==3
                                    Trans(3, 3) = 0;
                                    TransCost(3, 3) = 0;
                                else %m==4
                                    Trans(3, 4) = 0;
                                    TransCost(3, 4) = 0;
                                end
                                
                            else
                                if m==1
                                    Trans(4, 1) = Pw;
                                    TransCost(4, 1) = CW*Pw;
                                elseif m==2
                                    Trans(4, 2) = 0;
                                    TransCost(4, 2) = 0;
                                elseif m==3
                                    Trans(4, 3) = 0;
                                    TransCost(4, 3) = 0;
                                else %m==4
                                    Trans(4, 4) = (1-Pw);
                                    TransCost(4, 4) = CW*(1-Pw);
                                end
                                
                            end
                            
                            ECost = TransCost(n,m)*disc + Trans(n,m)*FutureOptNEC(m);
                            SumECost = SumECost + ECost;
                            
                        end % possible states at time t+1 (m)
                        SumECost = SumECost + (Xq + Xs + Xc)*disc;
                        
                        if SumECost < OptNEC
                            OptNEC = SumECost;
                            OptDec(n, 1, t, j) = Xq;
                            OptDec(n, 2, t, j) = Xs;
                            OptDec(n, 3, t, j) = Xc;
               
                        end
                        
                    end % allocations to control
                    
                end % allocations to surveillance
                
            end % allocations to quarantine
            
            FutureOptNEC(n) = OptNEC;
            OptCost(n,t,j) = OptNEC;
        end
        
    end
        
end
    

for k = 1:T
    AbsDec(1, k) = OptDec(1, 1, k, 1); %OptDec(state, action, time, discount rate)
    AbsDec(2, k) = OptDec(1, 2, k, 1);
    AbsDec(3, k) = OptDec(1, 3, k, 1);
    AbsDec(4, k) = OptDec(1, 1, k, 2);
    AbsDec(5, k) = OptDec(1, 2, k, 2);
    AbsDec(6, k) = OptDec(1, 3, k, 2);
end


AbsDec
time = 1:1:T;
plot(time, AbsDec,'LineWidth',2)
xlabel('Probability of spread (g)','FontSize',16)
ylabel('Investment (millions $)','FontSize',16)
legend('Quarantine, r = 0', 'Surveillance, r = 0', 'Removal, r = 0','Quarantine, r = 0.1', 'Surveillance, r = 0.1', 'Removal, r = 0.1')