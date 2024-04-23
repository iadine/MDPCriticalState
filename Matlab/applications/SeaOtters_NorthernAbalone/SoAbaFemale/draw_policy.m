% function draw_policy(Pi_Best,fname,liste)
% Action:
%   Graphical representation of a strategy (policy).
% Input:
%   Pi_Best: matrix representing the strategy,
%   fname:  string defining the FR and experts,
%   liste:  liste of action from function explore_Q
% Output:
%   Simple graph
%
% Author: iadine.chades@csiro.au

function draw_policy(Pi_Best,fname,liste)
global PARAM_MDP
global PARAM_ABALONE

unit=PARAM_ABALONE.discUnit;
nbs_so=PARAM_MDP.nbs_so;
nbs_aba=PARAM_MDP.nbs_aba;
x=0:unit:PARAM_ABALONE.k;
y=0:5:100;

Plan=zeros(nbs_so,nbs_aba);

for i=1:nbs_so
    for j=1:nbs_aba
        s=seeIndex([j-1 i-1]);
        Plan(i,j)=Pi_Best(s);
    end
end

A=[
    0.8 0.8 0.8
    1 0 0;
    0 1 0;
    0 0 0;
    0.5 0.5 0;
    ];
figure('color','white','Name',['Strategy computed for ',fname]);
image(x,y,Plan,'CDataMapping','scaled');
colormap(A);
caxis([0 5]);

hcb=colorbar('XTickLabel',...
    {' ','Nothing (N)',' ','Intro(I)',' '...
    'Anti-poach(A)','','Removal(R)',' ','1/2(A+R)'});
set(hcb,'YTickMode','manual','location','Southoutside')

ylabel 'Sea otter (%K)'
xlabel 'Abalone density (m^{-2})'
hold on
if isempty(liste) ~= 1
    for i=1:nbs_so
        for j=1:nbs_aba
            s=seeIndex([j-1 i-1]);
            if sum(liste(s,:))>1
                text((j-1)*unit,(i-1)*5,'*','color','white',...
                    'verticalalignment','Cap','horizontalalignment',...
                    'center','Fontweight','bold');
                hold on
            end
            
        end
    end
end
DirFig='Fig/';
saveas(gcf,[DirFig,fname], 'fig');
end

