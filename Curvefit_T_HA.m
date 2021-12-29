function [inf,sup,c_T,c_HA,t1,coefTHA,Degradat]=Curvefit_T_HA(u)
format long
% close all
c_T=0.1*u;
c_HA=1-c_T;
T=xlsread('T_CurveFit');
HA=xlsread('HA_CurveFit',3);
p=polyfit(T(:,1),T(:,2),3);
xT=T(:,1);
yT=T(:,2);
xHA=HA(:,1);
yHA=HA(:,2);
xTHA=c_T*xT+c_HA*xHA;
yTHA=c_T*yT+c_HA*yHA;

inf=4.05*c_T+6.45;
sup=5.5*c_T+12;
xyT=linspace(10.5,17.5,100);
xyHA=linspace(6.45,12,100);
xyTHA=linspace(inf,sup,100);
ft1 = fittype( 'exp1' );
ft2 = fittype( 'exp1' );
ft3 = fittype( 'exp1' );
% Fit model to data.
[fitresult1, gof1] = fit( xT, yT, ft1 );
[fitresult2, gof2] = fit( xHA, yHA, ft2 );
[fitresult3, gof3] = fit( xTHA, yTHA, ft3 );
% HApoly=polyval(p,HA(:,1));
% % f1=figure(1);
% % plot(xyT,fitresult1(xyT),'Color','k')
% % hold on
Degradat=[0.8*c_T 0.3*c_HA 0.6*c_HA 0.6];
% % plot(xyHA,fitresult2(xyHA),'Color','k')
% % plot(xyTHA,fitresult3(xyTHA),'Color',Degradat,'LineWidth',1.2)
% % coefT=coeffvalues(fitresult1);
% % coefHA=coeffvalues(fitresult2);
coefTHA=coeffvalues(fitresult3);
% % hold on
% % DarkRed1=[0.7 0 0];
% % DarkRed2=[0.7 0.2 0];
% % DarkRed3=[0.7 0.4 0];
% % DarkRed4=[0.7 0.6 0];
% % aT=scatter(T(1:10,1),T(1:10,2),'.','MarkerEdgeColor',DarkRed1);
% % bT=scatter(T(11:20,1),T(11:20,2),'.','MarkerEdgeColor',DarkRed2);
% % cT=scatter(T(21:30,1),T(21:30,2),'.','MarkerEdgeColor',DarkRed3);
% % dT=scatter(T(31:40,1),T(31:40,2),'.','MarkerEdgeColor',DarkRed4);
% % aT.SizeData=200;
% % bT.SizeData=200;
% % cT.SizeData=200;
% % dT.SizeData=200;
% % 
% % aHA=scatter(HA(1:10,1),HA(1:10,2),'.','MarkerEdgeColor',[0 0 0.7]);
% % bHA=scatter(HA(11:20,1),HA(11:20,2),'.','MarkerEdgeColor',[0 0.2 0.7]);
% % cHA=scatter(HA(21:30,1),HA(21:30,2),'.','MarkerEdgeColor',[0 0.4 0.7]);
% % dHA=scatter(HA(31:40,1),HA(31:40,2),'.','MarkerEdgeColor',[0 0.6 0.7]);
% % aHA.SizeData=200;
% % bHA.SizeData=200;
% % cHA.SizeData=200;
% % dHA.SizeData=200;
% % 
% % 
t1=fitresult3(xyTHA);
% % line([0 10.5],[3.8e-6 3.8e-6],'Color','k')
% % line([17.5 22],[1.92e-6 1.92e-6],'Color','k')
% % line([0 6.45],[4e-6 4e-6],'Color','k')
% % line([12 22],[1.859e-6 1.859e-6],'Color','k')
% % line([0 inf],[t1(1) t1(1)],'Color',Degradat,'LineWidth',1.2)
% % line([sup 22],[t1(end) t1(end)],'Color',Degradat,'LineWidth',1.2)
% % xticks(5:1:21)
% % ylim([1.5e-6 4.5e-6])
% % xlim([5 22])
% % xlabel('T/HA')
% % ylabel('α')
% % % xline(17.5,':')
% % % xline(10.5,':')
% % % xline(6.4,':')
% % % xline(12,':')
% % 
% % xline(inf,':')
% % xline(sup,':')
% % legend('','','Exponencial (T,HA)','','','','','','','','','','','','','')
% % 
% % text(16,4.1e-6,['Pes relatiu T=',num2str(100*c_T),'%'],'FontSize',8)
% % text(16,3.9e-6,['Pes relatiu HA=',num2str(100*c_HA),'%'],'FontSize',8)
% % hold off
% 
% print(f1,['CurvF_T_HA',num2str(u*5)],'-dpng','-r600')
% end
% close all
% video=VideoWriter('Gif_CurvF_T_HA.avi');
% open(video);
% for i=0:20
%    for v=1:30
%    I=imread(['CurvF_T_HA',num2str(i*5),'.png']);
%    writeVideo(video,I)
%    end
% end
% 
% close(video)


