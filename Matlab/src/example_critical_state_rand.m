function example_critical_state_rand()

% Define the parameters
S = 10;
A = 4;
discount = 0.96;
precision = 0.001;


% Generate random MDP example
[rand_P, rand_R] = mdp_example_rand(S, A);


% Define the path
c_dir =pwd();
pathToCheck = [c_dir, '\MDPToolbox\'];

% Check if the path is already added
if ~contains(userpath, pathToCheck)
    % If not, add the path
    addpath(pathToCheck);
end

% call function mdp_critical_state
[LP, LPplus,LPminus] = mdp_critical_state(rand_P, rand_R, discount, precision);

% Plot figure
figure('Color','white');
plot(1:S, LPminus*100, 'ko');
hold on;
plot(1:S, LPplus*100, 'ro');
for i = 1:S
    line([i, i], [LPminus(i)*100, LPplus(i)*100], 'Color', 'blue');
end
xlabel('States');
ylabel('Loss (vs optimal) (%)');
%title('Plot of LP- and LP+');
legend('LP- (best case)', 'LP+ (worst case)', 'Location', 'northeast');
hold off;
