close all
format long
HA=xlsread('HA_CurveFit',3);
p=polyfit(HA(:,1),HA(:,2),3);
x=HA(:,1);
y=HA(:,2);
xy=linspace(6.42,12,100);
ft = fittype( 'exp1' );
% Fit model to data.
[fitresult, gof] = fit( x, y, ft );
% HApoly=polyval(p,HA(:,1));
plot(xy,fitresult(xy),'Color','k')
coef=coeffvalues(fitresult);
% text(9.95,2.5e-6,['\leftarrow','alfa=',num2str(coef(1)),'·EXP(',num2str(coef(2)),'·HA)'],'FontSize',8)
% text(5.1,4.1e-6,'alfa=4.023e-6','FontSize',8)
% text(13,1.93e-6,'alfa=1.859e-6','FontSize',8)
hold on
DarkBlue1=[0 0 0.7];
DarkBlue2=[0 0.4 0.7];
DarkBlue3=[0 0.6 0.7];
DarkBlue4=[0 0.8 0.7];
a=scatter(HA(1:10,1),HA(1:10,2),'.','MarkerEdgeColor',DarkBlue1);
b=scatter(HA(11:20,1),HA(11:20,2),'.','MarkerEdgeColor',DarkBlue2);
c=scatter(HA(21:30,1),HA(21:30,2),'.','MarkerEdgeColor',DarkBlue3);
d=scatter(HA(31:40,1),HA(31:40,2),'.','MarkerEdgeColor',DarkBlue4);
a.SizeData=200;
b.SizeData=200;
c.SizeData=200;
d.SizeData=200;
line([5 6.42],[4.023e-6 4.023e-6],'Color','k')
line([12 15],[1.859e-6 1.859e-6],'Color','k')
xticks(5:1:15)
ylim([1.5e-6 4.5e-6])
xlabel('HA (g/m3)')
ylabel('α')
xline(12,':')
xline(6.4,':')
legend('','HA80','HA161','HA500','HA1500','','','','')
hold off




