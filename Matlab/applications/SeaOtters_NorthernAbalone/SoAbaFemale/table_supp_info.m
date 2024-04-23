function table_supp_info(filedate,EOA)
% for a given k and r provides table
% optimal
% antipoaching
% >60%K remove
% >60%K 1/2 remove and antipoaching
% Do nothing
% example: >> table_supp_info('24-Jun-2011','A')

global PARAM_SIG_FR PARAM_HYP_FR PARAM_LINEAR_FR DIR_results
% EOA= Experts only, Optimal only, All.
DIR_results='Results\';
if EOA=='E' || EOA==1
    EOA=1;
    M=dlmread([DIR_results,...
        filedate,'_results_Expert1.txt'],'',[0,0, 0, 4]);
    isAPoachEfficient=M(1); maxKAba = M(2);maxrAba=M(3);KSO = M(4); rSO=M(5);
elseif EOA=='O'|| EOA==2
    EOA=2;
M=dlmread([DIR_results,...
        filedate,'_results_VI.txt'],'',[0,0, 0, 4]);
    isAPoachEfficient=M(1); maxKAba = M(2);maxrAba=M(3);KSO = M(4); rSO=M(5);
elseif EOA=='A' || EOA==3
    EOA=3;
    M=dlmread([DIR_results,...
        filedate,'_results_Expert1.txt'],'',[0,0, 0, 4]);
    isAPoachEfficient=M(1); maxKAba = M(2);maxrAba=M(3);KSO = M(4); rSO=M(5);
end

% set defaults for optional inputs

load_param(isAPoachEfficient,maxKAba,maxrAba,KSO,rSO);

% A represents the map color for different expert/strategies
A=[
    0 1 0.2;
    0.1 0.1 0.1;
    0.5 0.5 0;
    0.5 0.5 0.5;
    0.0 0.5 1;
    ];

% Black and white figures.
Abw=[
    0 0 0;
    0.5 0.5 0.5;
    0 0 0;
    0.5 0.5 0.5;
    0.5 0.5 0.5;
    ];
% Initialise the parameters
v=length(PARAM_LINEAR_FR);      % v is the number of linear FR

nExpert=4+1;                    % 4 "experts" and an optimal cookie
res_exp=zeros(nExpert,v+4,8);   % v linear FR + 2 sig + 2 hyp

if EOA==1 || EOA==3
    % Retrieving the simulation results from files
    res_exp(2,:,:)=get_results([DIR_results,filedate,'_results_Expert1.txt']);
    res_exp(3,:,:)=get_results([DIR_results,filedate,'_results_Expert2.txt']);
    res_exp(4,:,:)=get_results([DIR_results,filedate,'_results_Expert3.txt']);
    res_exp(5,:,:)=get_results([DIR_results,filedate,'_results_Expert4.txt']);
end
if EOA==2 || EOA==3
    % Expert 5 is the optimal cookie
    res_exp(1,:,:)=get_results([DIR_results,filedate,'_results_VI.txt']);
end

% organise the res_exp matrices so that the data is in suitable format for
% the supp info table of Chades et al.

tab_s2=zeros(v+4*5,8);    % v+4 FRs; 5 experts; 8 fields
for i=1:v+4 % assumed 'A' input parameter
    for j=1:5 % experts
        tab_s2((i-1)*5+j,1:2)=res_exp(j,i,5:6);
        tab_s2((i-1)*5+j,3:4)=res_exp(j,i,3:4);
        tab_s2((i-1)*5+j,5:6)=res_exp(j,i,7:8);
        tab_s2((i-1)*5+j,7:8)=res_exp(j,i,1:2);
    end
      
end

dlmwrite('Results\tab_s2.txt', tab_s2);

