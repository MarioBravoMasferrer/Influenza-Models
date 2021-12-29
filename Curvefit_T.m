close all
format long
T=xlsread('T_CurveFit');
p=polyfit(T(:,1),T(:,2),3);
x=T(:,1);
y=T(:,2);
xy1=linspace(10.5,17.5,100);
xy2=linspace(9.5,19,100);
ft1 = fittype( 'exp1' );
ft2 = fittype( 'poly2' );
% Fit model to data.
[fitresult1, gof1] = fit( x, y, ft1 );
% HApoly=polyval(p,HA(:,1));
figure()
plot(xy1,fitresult1(xy1),'Color','k')
hold on
coef=coeffvalues(fitresult);
% text(9.95,2.5e-6,['\leftarrow','alfa=',num2str(coef(1)),'·EXP(',num2str(coef(2)),'·HA)'],'FontSize',8)
% text(5.1,4.1e-6,'alfa=4.02e-6','FontSize',8)
% text(13,1.93e-6,'alfa=1.9e-6','FontSize',8)
hold on
DarkRed1=[0.7 0 0];
DarkRed2=[0.7 0.2 0];
DarkRed3=[0.7 0.4 0];
DarkRed4=[0.7 0.6 0];
a=scatter(T(1:10,1),T(1:10,2),'.','MarkerEdgeColor',DarkRed1);
b=scatter(T(11:20,1),T(11:20,2),'.','MarkerEdgeColor',DarkRed2);
c=scatter(T(21:30,1),T(21:30,2),'.','MarkerEdgeColor',DarkRed3);
d=scatter(T(31:40,1),T(31:40,2),'.','MarkerEdgeColor',DarkRed4);
a.SizeData=200;
b.SizeData=200;
c.SizeData=200;
d.SizeData=200;
line([0 10.5],[3.8e-6 3.8e-6],'Color','k')
line([17.5 22],[1.9e-6 1.9e-6],'Color','k')
xticks(5:1:21)
ylim([1.5e-6 4.5e-6])
xlim([7 21])
xlabel('T (ºC)')
ylabel('α')
xline(17.5,':')
xline(10.5,':')
legend('','T80','T161','T500','T1500','','','','')
hold off




