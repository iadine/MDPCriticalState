%We have two plots, and we DON'T know how many golf tees are in each plot.
%Objective is to maximise the total number of captures

tic

T = 10; %Number of time steps

B = 2; %Budget - search time/effort to allocate between plots in each time step (in hours)
x1values = 0:1:B; %Vector of possible time allocated to plot 1 in each time setp

maxX1 = B*T;

maxN1 = 30; %Maximum number of plants that could be at site 1
maxN2 = 30; %ditto site 2

lambda1 = 0.07392; %Search efficiency in plot 1 
lambda2 = 0.03696; %Search efficiency in plot 2 
%Search efficiencies must be calibrated against the units used for the budget, currently hours per hectare.

buffer = 0.000000000001;

V = zeros(length(x1values),1); %value function for all decisions - expectated num captures in current timestep (for a particular state combination)
FutureVOpt = zeros(maxN1+1, maxN2+1, maxX1+1); %Value function for optimal decisions at time t+1
CurrentVOpt = zeros(maxN1+1, maxN2+1, maxX1+1); %Value function for optimal decisions at time t
PassiveDec = zeros(maxN1+1, maxN2+1, maxX1+1, T);
Values = zeros(maxN1+1, maxN2+1, maxX1+1, B+1);

posterior1_start = zeros(maxN1+1, maxN1+1, maxX1+1);
posterior1_now = zeros(maxN1+1, maxN1+1, maxX1+1);
posterior2_start = zeros(maxN2+1, maxN2+1, maxX1+1);
posterior2_now = zeros(maxN2+1, maxN2+1, maxX1+1);

binomial_prob1 = zeros(maxN1+1, maxN1+1, B+1);
binomial_prob2 = zeros(maxN2+1, maxN2+1, B+1);

prob1 = zeros(maxN1+1, maxX1+1, B+1, maxN1+1);
prob2 = zeros(maxN2+1, maxX1+1, B+1, maxN2+1);

%This prints out the optimal decisions, rewards, and capture probabilities in nicely formatted, readable files
xid = fopen('passive_dec_sim.dat', 'w'); %Optimal decisions in a format to be read into the simulation *Note that it read in backwards, from t = T:-1:1

yid = fopen('passive_optimal_dec.dat', 'w'); %Optimal decisions in a format to look at in Excel
fprintf(yid, 'Optimal decisions and their value\r\n\r\n');
fprintf(yid, '\r\n%s %d\r\n','Budget = ', B);
fprintf(yid, '\r\n%s %d, %s %d\r\n','maxN1 = ', maxN1, 'maxN2 = ', maxN2);
fprintf(yid, '\r\n%s %d, %s %d\r\n','lambda1 = ', lambda1, 'lambda2 = ', lambda2);

zid = fopen('passive_trans_probs.dat', 'w'); %

%Create matrices for the binomial probabilities of finding c plants out of n plants, for different values of c, n, and x.
%NOTE THAT FOR EFFICIENCY, THIS ASSUMES maxN1 = maxN2 AND maxX1 = maxX2
jid = fopen('passive_binomial_probs.dat', 'w');
fprintf(jid, 'Binomial probabilities\r\n\r\n');
fprintf(jid, '%s\r\n', 'Number in plot, Number found, Time searched, Prob plot 1, Prob plot 2');
for nindex = 1:maxN1+1 %Possible number of plants in a plot
    n = nindex-1;
    
    for cindex = 1:nindex %Possible number of plants found in a plot in a single time step
        c = cindex-1;
        
        for xindex = 1:B+1
            x = xindex-1;
        
            p1 = 1 - exp(-lambda1*x); %prob detecting an individual plant in site 1
            p2 = 1 - exp(-lambda2*x); %prob detecting an individual plant in site 2
        
            binomial_prob1(cindex, nindex, xindex) = calc_combinatorial(n,c).*(p1.^c).*((1-p1).^(n-c));
            binomial_prob2(cindex, nindex, xindex) = calc_combinatorial(n,c).*(p2.^c).*((1-p2).^(n-c));
            fprintf(jid, '%d, %d, %d, %f, %f\r\n', n, c, x, binomial_prob1(cindex, nindex, xindex), binomial_prob2(cindex, nindex, xindex));
        end
    end
end

%Create matrices of posterior probability distributions for the number originally in the site and the number remaining at time t
%NOTE THAT FOR EFFICIENCY, THIS ASSUMES maxN1 = maxN2 AND maxX1 = maxX2
for Cindex = 1:maxN1+1 %Possible number captured in plot 1 across entire management period
    C = Cindex-1; %Number already found and removed
    
    for Xindex = 1:maxX1+1 %Possible cumulative time spent searching plot 1 so far
        X = Xindex-1; %Time already searched

        %Calculating posterior prob distribution for the number of plants initially in site 1 (N0)
        for i1index = Cindex:maxN1+1 %Must have been at least C plants (number already captured)
            i1 = i1index - 1; %possible numbers of plants originally in site 1
            posterior1_start(i1index, Cindex, Xindex) = calc_combinatorial(i1, C).*exp(-lambda1.*X.*i1);
        end
        Sum1 = sum(posterior1_start(:,Cindex,Xindex));
        %Normalising so probabilities sum to one
        for j1index = Cindex:maxN1+1
            j1 = j1index - 1;
            %This is the probability that there were j1 plants at the start, given C and X
            %It's a distribution over the possible values of j1
            posterior1_start(j1index, Cindex, Xindex) = posterior1_start(j1index, Cindex, Xindex)./Sum1;
            
            %Now this is the probability that there are num1 plants NOW, given C and X
            %It's the previous distribution, shifted down by C
            num1 = j1 - C; %possible numbers of plants in site 1 now
            num1index = num1+1;
            posterior1_now(num1index,Cindex, Xindex) = posterior1_start(j1index,Cindex, Xindex);
        end
        
        if lambda2 == lambda1
            posterior2_start(:, Cindex, Xindex) = posterior1_start(:, Cindex, Xindex);
            posterior2_now(:, Cindex, Xindex) = posterior1_now(:, Cindex, Xindex);           
        else
            
            %Do the same for site 2
            for i2index = Cindex:maxN2+1 %Must have been at least C plants (number already captured)
                i2 = i2index - 1; %possible numbers of plants originally in site 2
                posterior2_start(i2index, Cindex, Xindex) = calc_combinatorial(i2, C).*exp(-lambda2.*X.*i2);
            end
            Sum2 = sum(posterior2_start(:,Cindex,Xindex));
            %Normalising so probabilities sum to one
            for j2index = Cindex:maxN2+1
                j2 = j2index - 1;
                %This is the probability that there were j2 plants at the start, given C and X
                %It's a distribution over the possible values of j2
                posterior2_start(j2index, Cindex, Xindex) = posterior2_start(j2index,Cindex, Xindex)./Sum2;
                
                %Now this is the probability that there are num2 plants NOW, given C and X
                %It's the previous distribution, shifted down by C
                num2 = j2 - C; %possible numbers of plants in site 2 now
                num2index = num2+1;
                posterior2_now(num2index,Cindex, Xindex) = posterior2_start(j2index,Cindex, Xindex);
            end
        end  
    end
end                

%This next section prints the posterior probability distributions to a file
kid = fopen('passive_posteriors.dat', 'w');
fprintf(kid, 'Posterior probability distributions\r\n\r\n');
fprintf(kid, '%s\r\n', 'C1, X1, Number, Prob1 Start, Prob1 Now');
for Cindex = 1:maxN1+1 %Possible number captured in plot 1 across entire management period
    C = Cindex - 1;
    
    for Xindex = 1:maxX1+1 %Possible cumulative time spent searching plot 1 so far
        X = Xindex - 1;
        
        for iindex = 1:maxN1+1
            i = iindex - 1;
            
            fprintf(kid, '%d, %d, %d, %f, %f\r\n', C, X, i, posterior1_start(iindex, Cindex, Xindex), posterior1_now(iindex, Cindex, Xindex));
            
        end 
        
    end
end
fprintf(kid, '\r\n\r\n');   
fprintf(kid, '%s\r\n', 'C2, X2, Number, Prob2 Start, Prob2 Now');
for Cindex = 1:maxN2+1 %Possible number captured in plot 1 across entire management period
    C = Cindex - 1;
    
    for Xindex = 1:maxX1+1 %Possible cumulative time spent searching plot 1 so far
        X = Xindex - 1;
        
        for iindex = 1:maxN2+1
            i = iindex - 1;
            
            fprintf(kid, '%d, %d, %d, %f, %f\r\n', C, X, i, posterior2_start(iindex, Cindex, Xindex), posterior2_now(iindex, Cindex, Xindex));
            
        end 
        
    end
end


%Create matrix of final rewards 
for C1index = 1:maxN1+1 %Possible number captured in plot 1 across entire management period
    C1 = C1index-1; %C1,t-1
    
    for C2index = 1:maxN2+1 %Possible number captured in plot 2 across entire management period
        C2 = C2index-1; %C2,t-1
        
        for X1index = 1:maxX1+1 %Possible cumulative time spent searching plot 1 so far
            X1 = X1index-1; %X1
            X2 = maxX1 - X1; %X2
            
            if X1 == 0 && C1 > 0 %If we haven't looked in site 1, there can't have been captures there
                FutureVOpt(C1index, C2index, X1index) = FutureVOpt(C1index, C2index, X1index)*NaN; %Set as Not a Number
                
            elseif X2 == 0 && C2 > 0 %If we haven't looked in site 2, there can't have been captures there
                FutureVOpt(C1index, C2index, X1index) = FutureVOpt(C1index, C2index, X1index)*NaN; %Set as Not a Number
                    
            else
                FutureVOpt(C1index, C2index, X1index) = C1+C2; %Final reward is total number of captures    
            
            end
        end
    end
end


%SDP - Loop over all time steps
for t = T:-1:1
    t
    tic
    maxsearch = B*(t-1); %The maximum amount of total search time so far (up until start of time t)
    
    fprintf(yid, '\r\n%s %d\r\n','time = ', t);
    fprintf(yid, '%s\r\n','X1, X2, C1, C2, V for decision 0, V for decision 1, V for decision 2, OptV, Optx1,');
    
    for C1index = 1:maxN1+1 %Possible number captured in plot 1 so far
        C1 = C1index-1; %C1,t-1
        
        for C2index = 1:maxN2+1 %Possible number captured in plot 2 so far
            C2 = C2index-1; %C2,t-1
            
            for X1index = 1:maxX1+1 %Possible cumulative time spent searching plot 1 so far
                X1 = X1index-1; %X1,t-1
                X2 = maxsearch - X1; %X2,t-1
                X2index = X2+1;
                
                if X1 > maxsearch %If X1 > max amount of search time so far, then state not possible
                    CurrentVOpt(C1index, C2index, X1index) = CurrentVOpt(C1index, C2index, X1index)*NaN; %Set as Not a Number
                    PassiveDec(C1index, C2index, X1index, t) = PassiveDec(C1index, C2index, X1index, t)*NaN;
                    
                elseif X1 == 0 && C1 > 0 %If we haven't looked in site 1, there can't have been captures there
                    CurrentVOpt(C1index, C2index, X1index) = CurrentVOpt(C1index, C2index, X1index)*NaN; %Set as Not a Number
                    PassiveDec(C1index, C2index, X1index, t) = PassiveDec(C1index, C2index, X1index, t)*NaN;
                
                elseif X2 == 0 && C2 > 0 %If we haven't looked in site 2, there can't have been captures there
                    CurrentVOpt(C1index, C2index, X1index) = CurrentVOpt(C1index, C2index, X1index)*NaN; %Set as Not a Number
                    PassiveDec(C1index, C2index, X1index, t) = PassiveDec(C1index, C2index, X1index, t)*NaN;
                
                else
                    %If none of the previous conditions apply, find optimal decision
                    for x1index = 1:length(x1values) %Possible search time allocated to plot 1 in THIS timestep
                        
                        x1 = x1values(x1index); %x1,t
                        x2 = B - x1; %What's left over is allocated to plot 2 (x2,t)
                        x2index = x2 + 1;
                        
                        newX1 = X1 + x1; %Update the cumulative amount of search time
                        newX2 = X2 + x2;
                                              
                        Sum_expcap1 = 0;
                        for c1index = 1:(maxN1-C1+1) %Possible number of captures in site 1 in THIS time step 
                            c1 = c1index-1; %c1,t
                            prob1(c1index,1) = binomial_prob1(c1index, 1:(maxN1-C1+1), x1index)*posterior1_now(1:(maxN1-C1+1), C1index, X1index);
                            Sum_expcap1 = Sum_expcap1 + c1.*prob1(c1index,1);
                        end
                        Expcap1 = Sum_expcap1;
                                
                        Sum_expcap2 = 0;
                        for c2index = 1:(maxN2-C2+1) %Possible number of captures in site 2 in THIS time step
                            c2 = c2index-1; %c2,t
                            prob2(c2index,1) = binomial_prob2(c2index, 1:(maxN2-C2+1), x2index)*posterior2_now(1:(maxN2-C2+1), C2index, X2index);
                            Sum_expcap2 = Sum_expcap2 + c2.*prob2(c2index,1);
                        end
                        Expcap2 = Sum_expcap2;
                                
                        V(x1index,1) = Expcap1 + Expcap2;

                    end
                    
                    Diff01 = abs(V(1) - V(2));
                    Diff02 = abs(V(1) - V(3));
                    Diff12 = abs(V(2)- V(3));
                    
                    %This function finds and stores the value and index of the maximum value in V
                    [CurrentVOpt(C1index, C2index, X1index), XOptIndex] = max(V);
                    XOpt = x1values(XOptIndex); %Finds the decision with that index
                    
                    %Flag if value functions for different decisions are really similar
                    if XOpt == 0
                        if Diff01 < buffer
                            if Diff02 < buffer
                                XOpt = 7; %120, all decisions equal
                            else
                                XOpt = 5; %10 Decisions 0 and 1 are equal
                            end
                        elseif Diff02 < buffer
                            XOpt = 4; %20 Decisions 0 and 2 are equal
                        end
                    end
                    if XOpt == 1
                        if Diff01 < buffer
                            if Diff12 < buffer
                                XOpt = 7; %120, all decisions equal
                            else
                                XOpt = 5; %10 Decisions 0 and 1 are equal
                            end
                        elseif Diff12 < buffer
                            XOpt = 6; %12 Decisions 1 and 2 are equal
                        end
                    end
                    if XOpt == 2
                        if Diff12 < buffer
                            if Diff02 < buffer
                                XOpt = 7; %120, all decisions equal
                            else
                                XOpt = 6; %12 Decisions 1 and 2 are equal
                            end
                        elseif Diff02 < buffer
                            XOpt = 4; %20 Decisions 0 and 2 are equal
                        end
                    end
                    %
                    PassiveDec(C1index, C2index, X1index, t) = XOpt;
                    Values(C1index, C2index, X1index, :) = V;
                    
                    fprintf(yid, '%d, %d, %d, %d, %f, %f, %f, %f, %d\r\n', X1, X2, C1, C2, V(1), V(2), V(3), CurrentVOpt(C1index, C2index, X1index), XOpt);
                    V = zeros(length(x1values),1); %clear matrix                    
                end   
            end            
        end        
    end

    %print decisions to file in format to be read in to simulation
    for iindex = 1:maxN1+1 %C1index
        
        for jindex = 1:maxN2+1 %C2index
            
            for kindex = 1:maxsearch+1 %X1index
                
                fprintf(xid, '%d, ', PassiveDec(iindex, jindex, kindex, t));
            end
            fprintf(xid, '\r\n');    
        end
    end
    
    FutureVOpt = CurrentVOpt;
    toc
end

Cvalues = 0:30; %corrects axes
createfigure(Cvalues, PassiveDec(:, :, 2, 2)); 

