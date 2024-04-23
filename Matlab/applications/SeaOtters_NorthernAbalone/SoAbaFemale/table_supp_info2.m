function table_supp_info2(filedate,EOA)
% provides table s3 - use results_analysis_SA
% A faire
% optimal
% antipoaching
% >60%K remove
% >60%K 1/2 remove and antipoaching
% Do nothing
% example: >> table_supp_info2('19-Aug-2011','O')

global PARAM_SIG_FR PARAM_HYP_FR PARAM_LINEAR_FR DIR_results
% EOA= Experts only, Optimal only, All.
DIR_results='Results\';
if EOA=='O' || EOA==2
    EOA=2;
    M=dlmread([DIR_results,...
        filedate,'_results_SA_VI.txt'],'',[0,0, 0, 4]);
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
new_k=[3.34-0.2,3.34 , 3.34+0.2];
new_r=[1.6-0.05, 1.6, 1.6+0.05];

nSA=length(new_k)*length(new_r);
res_SA=zeros(nSA,v,8);   % v linear FR + 2 sig + 2 hyp
X=get_results([DIR_results,filedate,'_results_SA_VI.txt']);
for i=1:nSA
    for ifr=1:v
        res_SA(i,ifr,:)=X(v*(i-1)+ifr,:);
    end
end

% organise the res_SA matrices so that the data is in suitable format for
% the supp info table of Chades et al.

tab_s3=zeros(nSA*v,8);    % v FRs; 8 fields

for i=1:v % assumed 'A' input parameter
    for j=1:nSA % parameters
        tab_s3((i-1)*nSA+j,1:2)=res_SA(j,i,5:6);
        tab_s3((i-1)*nSA+j,3:4)=res_SA(j,i,3:4);
        tab_s3((i-1)*nSA+j,5:6)=res_SA(j,i,7:8);
        tab_s3((i-1)*nSA+j,7:8)=res_SA(j,i,1:2);
    end

end

dlmwrite('Results\tab_s3.txt', tab_s3);

