% function All=get_results(name)
% Action:
%   Reads and extracts results for a given file name that identifies a
%   strategy and its performance. Called when results are analysed.
%   
% Input:
%   name: filename
% Output:
%   All: [score, std_sc,aba,std_aba,aba100,std_aba100,so,std_so];
% Side effect:
%   The structure of the file is highly dependent on the 
%   evaluate_strategy_mean function.
% Author:
%   iadine.chades@csiro.au

function All=get_results(name)

X=dlmread(name,'',3,0);     % file dependent
score=X(:,1);
std_sc=X(:,2);
aba=X(:,5);
std_aba=X(:,6);
so=X(:,3);
std_so=X(:,4);
aba100=X(:,7);
std_aba100=X(:,8);
All=[score, std_sc,aba,std_aba,aba100,std_aba100,so,std_so];
end

