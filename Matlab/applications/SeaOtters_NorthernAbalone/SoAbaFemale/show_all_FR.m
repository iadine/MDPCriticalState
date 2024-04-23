function show_all_FR()
global PARAM_LINEAR_FR

% Initialising/loading parameters
disp('-> Loading parameters');
load_param(0,3.34,1.6,4073,0.191);
disp('<- Parameters loaded');

figure('color','white');

v=length(PARAM_LINEAR_FR);
for iFR=1:2
    for iFRnum=1:2
        show_FR(iFR,iFRnum) % Plot FR
    end
end
for iFRnum=1:v
    show_FR(3,iFRnum) % Plot FR
end