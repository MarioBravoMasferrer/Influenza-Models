% function [error]=ALFA_T_14(x0)
%% Optim model alfa f(HA)
close all

%% ALFA F(HA)
%% Paràmetres
v_anys = [2010,2011,2012,2013,2014,2015,2016,2017,2018,2019];
[num]=xlsread('CasosR15_18.xlsx');                     
T=xlsread('T_BCN');        
T=T(1:end,:);
er=[];
DeltaT = 1;           %dia
Tini=1;               %dies
Tfin=357;             %dies
beta = 1/4;           %1/dia
gamma = 1/7;          %1/dia
N=5e6;            %ind   

R2=[];p_val=[];R2G=[];
% F=[0.01665;0.01612;0.01702;0.0169;0.01782;0.0172;0.01824;0.01946;0.01786;0.01741]; % Error per sota e8
F=round([0.0166125841068557;0.0162029363068831;0.0188291383125575;0.0175934687370616;0.0181293096651430;0.0175829707652919;0.0190183461981882;0.0199382110038480;0.0198955736778908;0.0183563232925347],6);  %Error 2.6e7
% F=[0.0164751837986060;0.0159094987783852;0.0183723926459344;0.0175620213713224;0.0178739841448648;0.0174582324347061;0.0181918430505514;0.0200171516157985;0.0197000135194311;0.0180757301817864];
% F=[0.01215;0.01062;0.01202;0.0129;0.01302;0.010;0.01224;0.01546;0.01186;0.01241]; %F temperatura guai

% Io=[60;54;50;43;31;26;43;25;45;38]; % Error per sota e8
Io=round([78.9319564747923;42.7022154546790;2.73624751839554;20.8802239678661;10.9194952092757;11.5157984399925;35.5020031318091;8.28367040312476;1.99459359678951;7.98970302045288]); % Error 2.6e7

% F=round([0.0163229665760769;0.0159094271958896;0.0184654935404624;0.0175619954220037;0.0178739598967948;0.0174191675159815;0.0181916989234967;0.0200081834057623;0.0197000101623078;0.0180757082803716],6); % F temperatura aburrit
% Io=round([75.6816250949958;52.6795475524000;3.77060958179829;15.6750839530685;11.3649942892719;9.91276211594153;82.0529620374410;5.97743008315192;2.25714083474767;9.52605101289737]); % Io Temperatura aburrit
% Io=[60.1688059217661;51.8650528767824;4.22104838735055;15.4373351864635;11.1816964946058;9.21296990236278;80.7047463527931;5.81177579493136;2.21795994811254;9.37406398337219];
id=[];
id2=datetime();
Temps=Tini:DeltaT:Tfin;
Npassos=length(Temps);
TempsCR=Tini:7:Tfin;
Error=zeros(1:length(v_anys));
error=0;

R2T=[];

llindT=[10.5,17.5];
% llindT=[9.73,19];



for i = 1:length(v_anys)    
Any = v_anys(i);
S=zeros(Npassos,1);
E=zeros(Npassos,1); 
I=zeros(Npassos,1);
R=zeros(Npassos,1);
alfa=zeros(Npassos,1);
Y=Any-2009;
f=F(Y);
% Casos Reals
CasosReals=num(1:51,Y+1);


% Valors Inicials
S(1)=f*N;
E(1)=0;
I(1)=Io(Y);
R(1)=0;
alfa(1)=1.7e-06;



for t=2:Npassos

%% Alfa en funció de curvefitting
%% Coeficients T
k=t;
% p=[1.46177589893955e-08,-6.49948294228726e-07,8.97852950346182e-06]; 
p=[1.06231501502854e-05,-0.0979123978957006];
    if T(k-1,Y)<=llindT(2) && T(k-1,Y)>=llindT(1)
        m1= p(1);
        n1=p(2);
        alfa(k)=m1*exp(T(k-1,Y)*n1);
%         alfa(k)=m1*T(k-1,Y)^n1; 
%         alfa(k)=p(1)*T(k-1,Y)^2+p(2)*T(k-1,Y)+p(3);
        
    end
    
    if T(k-1,Y)<llindT(1)
    alfa(k)=3.8e-6;
    end
    if T(k-1,Y)>llindT(2)
    alfa(k)=1.9e-6; 
    end


    S(t)=S(t-1)-(alfa(t-1)*S(t-1)*I(t-1))*DeltaT;
    E(t)=E(t-1)+(alfa(t-1)*S(t-1)*I(t-1)-beta*E(t-1))*DeltaT;
    I(t)=I(t-1)+(beta*E(t-1)-gamma*I(t-1))*DeltaT;
    R(t)=R(t-1)+(gamma*I(t-1))*DeltaT;   


end
% [r,p]=corrcoef(I(1:7:end),CasosReals(1:51));

[~, idx]=max(I);
x=sort(I(1:7:floor(length(Temps(1:idx)))-6));
y=sort(CasosReals(1:idx/7));
mdl = fitlm(x,y);
r2=mdl.Rsquared.Ordinary;
R2(i)=mdl.Rsquared.Ordinary;
p_val(i)=mdl.Coefficients.pValue(2);

%% Gràfiques 
d1=cell(52,1);
for p=23:52
    d1(p-22)={num2str(p)};
end
for p=1:22
    d1(p+30)={num2str(p)};
end
f1 = figure(1);
f1.Position = [100 10 1100 1200];
% Gràfiques
Ini23=[7,6,4,3,2,1,6,5,4,3];
Red = 1/255*[220,50,0];
v=i;
if i>8
    v=i+1;
end
% d1 = datetime(Any,6,Ini23(Y)) + caldays(1:length(I));
subplot(3,4,v)
hold on
Infectats=plot(Temps,I(:),'k','LineWidth',1.2);
hold on
plot(TempsCR,CasosReals,'.','Color',Red,'MarkerSize',8)
set(gca,'xtick',1:28:length(I));
set(gca,'xticklabel',d1(1:4:end))
title([num2str(Any),'-',num2str(Any+1)])
xlabel("Setmanes de l'any");
if i==1 || i==5 || i==9
ylabel('Individus Infectats')
end
xlim([0 364])
ylim([0 17000])
set(gca,'ytick',0:2000:17000);
text(10,14500,['R2= ',num2str(R2(i))],'FontSize',7)
text(10,16000,['p-value= ',num2str(p_val(i))],'FontSize',7)
if i==1
legend('Model','Casos Reals','location','northwest');
end
hold off

if i==1
figure()
colororder({'k','k'})
plot(Temps,T(14:length(Temps)+13,Y),'Color',[1 0 0.3 0.9],'LineWidth',1.2)
xlabel("Setmanes de l'any")
ylabel("Temperatura (ºC)")
set(gca,'xtick',1:7:length(I));
set(gca,'xticklabel',d1(1:1:end))
hold on
yyaxis right
xlim([0 360])
plot(Temps,alfa(:),'Color','k','LineWidth',1.1)
ylabel('alfa')
legend('T (ºC)','Alfa')
hold off
end
%% ERROR
% o=0;
% for l=1:length(I)
%    if I(l)>200 && o==0
%        idx=l;
%        o=1;
%    end
% end
% id(i)=I(idx);
% id2(i)=d1(idx);


% o=0;
% for l=1:51
%    if CasosReals(l)>1500 && o==0
%        km=l;
%        o=1;
%    end
% end
% id2(i)=d1(km*7-6);
[~, idx]=max(I);
Mid=floor(length(Temps(1:idx))/7);
new_I_matrix=zeros(Mid,length(v_anys));


len=floor(length(Temps(1:idx))/7);
new_I_matrix(1:len,i)=I(1:7:floor(length(Temps(1:idx)))-6);

for j=1:Mid
    error = error + (new_I_matrix(j,i)-CasosReals(j))^2;
end
Error(i)=error;


%% CoefCorrel
figure(2)
subplot(3,4,i)
hold on
plot(x,y,'.','Color',[0.7 0 0],'MarkerSize',10)
plot(x,mdl.Fitted)
xlabel('Model')
ylabel('Casos')
text(200,9000,['R2=',num2str(mdl.Rsquared.Ordinary)],'FontSize',7)
ylim([0 10000])
hold off
R2G(i)=mdl.Rsquared.Ordinary;
end
 
%% Curve fitting alfas

error=sum(Error(:));
% T=10.5:0.1:17.5;
% func=[];
% for i=1:length(T)
% func(i)=m1*T(i)^n1;
% end
% figure()
% plot(10.5:0.1:17.5,func)