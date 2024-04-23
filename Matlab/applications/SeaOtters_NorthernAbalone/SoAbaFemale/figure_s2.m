function figure_s2()

global PARAM_ABALONE
load_param(0,3.34,1.6,4073,0.191)

rmax=derive_rmax(1.6);      % rmax=[1.05,1.2,1.4,1.6];  
kmax= derive_kabah(3.34); % k_aba_h=[0.837,1.67,2.5,3.34];  % new K

y=[rmax];
z=[kmax];
figure('color','white')
plot(0:3,y,'*k-')
hold on
yb=rmax+0.05;
plot(0:3,yb,'ok--')
yb=rmax-0.05;
plot(0:3,yb,'ok--')
set(gca,'Xtick',[0:3])
set(gca,'XTicklabel',{'Poor';'Medium';'Good';'Very good'});
Xlabel('Habitat quality');
Ylabel('Intrinsic growth rate');
legend('r assumed in (Chades et al, this paper)','r +/- 0.05 sensitivity analysis');
legend('boxoff');
box off

figure('color','white')

plot(0:3,z,'k*-')
hold on
zb=z+0.2;
plot(0:3,zb,'ko--')
hold on
zb=z-0.2;
plot(0:3,zb,'ko--')

set(gca,'Xtick',[0:3])
set(gca,'XTicklabel',{'Poor';'Medium';'Good';'Very good'});
Xlabel('Habitat quality');
Ylabel('Carrying capacity');
legend('K assumed in (Chades et al, this paper)','K +/- 0.2 sensitivity analysis');
legend('boxoff');
box off
end

function rmax=derive_rmax(maxrAba)  % default value is rmax~[1.05,1.2,1.4,1.6];  
 if maxrAba<1.05
     error('maxrAba has to be > 1.05')
 end
b=1.05;             % growth rate in poor habitat
a= (maxrAba-b)/3;   % assumed a linear relationship between different habitat type
rmax=[b,a+b,2*a+b,maxrAba];
end

function k_aba_h= derive_kabah(maxKAba) % k_aba_h=[0.837,1.67,2.5,3.34];  % new K
a=maxKAba/4;  % assumed a linear relationship between different habitat type
k_aba_h=[a, 2*a, 3*a, maxKAba];
end