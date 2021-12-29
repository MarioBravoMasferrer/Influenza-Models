% function [error]=Alfa_T_Girona(x0,u)
%% Optim model alfa f(HA)
close all

%% ALFA F(HA)
%% Paràmetres
v_anys = [2011,2012,2013,2014,2015,2016,2017,2018,2019];
% v_anys=v_anys(u);
[num]=xlsread('Grip Girona.xlsx');                     
T=xlsread('T_Girona');        
T=T(1:end,:);
er=[];
DeltaT = 1;           %dia
Tini=1;               %dies
Tfin=357;             %dies
beta = 1/4;           %1/dia
gamma = 1/7;          %1/dia
N=888467;
NB=5e6;            %ind   

R2=[];p_val=[];
% Io=num(1,3:end);
F=[0.0154034191814021;0.0164514722132001;0.0166080165575007;0.0174189463555523;0.0168441499241323;0.0174514079631569;0.0189879623758098;0.0187306488368297;0.0148018271755170];
Io=([4.09898111616654;1.9358812988172;1.92837795419100;1.23973823655097;0.234223497575379;3.67484433039413;0.758331521379005;0.242218946053037;30.4222547896818]);
% F=x0(1);
% Io=x0(2);
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
Y=Any-2010;
% Y=u;
f=F(Y);
% f=F;
% Casos Reals
CasosReals=num(1:51,Y+2);


% Valors Inicials
S(1)=f*N;
E(1)=0;
I(1)=Io(Y);
% I(1)=Io;
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
    
    alfa(k)=alfa(k)*NB/N;
    
    S(t)=S(t-1)-(alfa(t-1)*S(t-1)*I(t-1))*DeltaT;
    E(t)=E(t-1)+(alfa(t-1)*S(t-1)*I(t-1)-beta*E(t-1))*DeltaT;
    I(t)=I(t-1)+(beta*E(t-1)-gamma*I(t-1))*DeltaT;
    R(t)=R(t-1)+(gamma*I(t-1))*DeltaT;   


end
[r,p]=corrcoef(I(1:7:end),CasosReals(1:51));

R2(i)=r(1,2);
p_val(i)=p(1,2);

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
% % Gràfiques
Ini23=[7,6,4,3,2,1,6,5,4,3];
Red = 1/255*[220,50,0];
v=i;
if i>8
    v=i+1;
end

subplot(3,3,i)
hold on
Infectats=plot(Temps,I(:),'k','LineWidth',1.2);
hold on
plot(TempsCR,CasosReals,'.','Color',Red,'MarkerSize',8)
set(gca,'xtick',1:28:length(I));
set(gca,'xticklabel',d1(1:4:end))
title([num2str(Any),'-',num2str(Any+1)])
if i>6
xlabel("Setmanes de l'any");
end
if i==1 || i==4 || i==7
ylabel('Individus Infectats')
end
xlim([0 364])
ylim([0 2800])
set(gca,'ytick',0:500:3000);
text(10,2650,['R2= ',num2str(R2(i))],'FontSize',7)
text(10,2400,['p-value= ',num2str(p_val(i))],'FontSize',7)
if i==1
legend('Model','Casos Reals','location','northwest');
end
hold off

end
%% ERROR

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

R_2T=corrcoef(T(1:7:length(CasosReals)*7,Y),CasosReals);
R2T(i)=R_2T(1,2);



 

%% Optimització
% close all
% OUT=zeros(9,2);
% for u=1:9
% F=[0.0157263610600156;0.0171844976811344;0.0178235218329173;0.0181324729621484;0.0166419299472286;0.0173300240393561;0.0196197758219994;0.0182785448906585;0.0171558537780388];
% Io=[10,11,6,3,2,3,4,5,8];
% x0=[F(u),Io(u)];
% fun=@(x0)Alfa_T_Girona(x0,u);
% options = optimset('MaxFunEvals',100000,'MaxIter',100000);  % Optimització dels paràmetres de la funció 
% OUT(u,:)=fminsearch(fun,x0,options);
% end