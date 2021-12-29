% function [error]=Alfa_HA_Lleida(x0,u)
format long
close all
%% ALFA(HA)
%% Paràmetres
v_anys = [2010,2011,2012,2013,2014,2015,2016,2017,2018,2019];
% v_anys=v_anys(u);
[num]=xlsread('DR_Lleida.xlsx');                     
HA=xlsread('HA_Lleida');        
HA=HA(1:end,1:10);

DeltaT = 1;           %dia
Tini=1;               %dies
Tfin=357;             %dies
beta = 1/4;           %1/dia
gamma = 1/7;          %1/dia
N=365602;            %ind   
NB=5e6;
alfas=[];
F=[0.0166978737159249;0.0162725024403515;0.0173100410588718;0.0164288588908005;0.0218360043323922;0.0159045177502600;0.0174833193706885;0.0200497563567780;0.0191639895073391;0.0176856577053125];
Io=([0.424333037262160;0.259582994404641;0.107345559963256;1.03577355548904;0.00119842193847516;0.197130672433497;0.789585454491274;0.0238991390741969;0.127494101071115;0.103611341717825]);
% F=[0.0157905682864913,0.0153164471970013,0.0179776884551847,0.0168353145218735,0.0170066415007447,0.0161399935283330,0.0183365014233892,0.0185843401086952,0.0193544254824303,0.0179297124222972];  % Error 1.8e-8
% Io=round([86.5677108089879,76.1353112729644,2.78364165584862,16.6338446284347,46.8231764506672,40.6135955050504,39.7313976521193,16.4180257271306,3.91231285302538,5.39172980254877]*N/NB);
% Io=[6,6,1,1,3,3,3,1,1,1];
% F=x0(1);
% Io=x0(2);
Temps=Tini:DeltaT:Tfin;
Npassos=length(Temps);
TempsCR=Tini:7:Tfin;
Error=zeros(1:length(v_anys));
error=0;
llindHA=[11.9478332519531,6.42500152587891];


HAvals=zeros(51,21);

for i = 1:length(v_anys)    
Any = v_anys(i);
S=zeros(Npassos,1);
E=zeros(Npassos,1);
I=zeros(Npassos,1);
R=zeros(Npassos,1);
alfa=zeros(Npassos,1);
Y=Any-2009;
% Y=u;
f=F(Y);
% f=F;
% Casos Reals
CasosReals=num(1:51,Y);


% Valors Inicials
S(1)=f*N;
E(1)=0;
I(1)=Io(Y);
% I(1)=Io;
R(1)=0;
alfa(1)=1.9e-06;


for t=2:Npassos

%% Alfa en funció de curvefitting

k=t;
    if HA(k-1,Y)<=llindHA(1) && HA(k-1,Y)>=llindHA(2)
        p=[9.781e-6,-0.1383];
        alfa(k)=p(1)*exp(p(2)*HA(k-1,Y));

    end
    
    if HA(k-1,Y)>llindHA(1)
    alfa(k)=1.859e-06;

    end
    if HA(k-1,Y)<llindHA(2)
    alfa(k)=4.023e-06;

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
% Gràfiques
Ini23=[7,6,4,3,2,1,6,5,4,3];
Blue = 1/255*[0,150,200];
v=i;
if i>8
    v=i+1;
end
% d1 = datetime(Any,6,Ini23(Y)) + caldays(1:length(I));
subplot(3,4,v)
hold on
Infectats=plot(Temps,I(:),'k','LineWidth',1.2);
hold on
plot(TempsCR,CasosReals,'.','Color',Blue,'MarkerSize',8)
set(gca,'xtick',1:28:length(I));
set(gca,'xticklabel',d1(1:4:end))
title([num2str(Any),'-',num2str(Any+1)])
xlabel("Setmanes de l'any");
if i==1 || i==5 || i==9
ylabel('Individus Infectats')
end
xlim([0 364])
ylim([0 1500])
set(gca,'ytick',0:250:1500);
text(10,13500,['R2= ',num2str(R2(i))],'FontSize',7)
text(10,15000,['p-value= ',num2str(p_val(i))],'FontSize',7)
if i==1
legend('Model','Casos Reals','location','northwest');
end
hold off

[~, idx] = max(I);


%% ERROR

Mid=floor(length(Temps(1:idx))/7);
new_I_matrix=zeros(Mid,length(v_anys));


len=floor(length(Temps(1:idx))/7);
new_I_matrix(1:len,i)=I(1:7:floor(length(Temps(1:idx)))-6);

for j=1:Mid
    error = error + (new_I_matrix(j,i)-CasosReals(j))^2;
end
Error(i)=error;



end


error=sum(Error(:));

%% Optimització
% close all
% OUT=zeros(10,2);
% for u=1:10
% F=[0.0157905682864913,0.0153164471970013,0.0179776884551847,0.0168353145218735,0.0170066415007447,0.0161399935283330,0.0183365014233892,0.0185843401086952,0.0193544254824303,0.0179297124222972];
% Io=[6,6,1,1,3,3,3,1,1,1];
% x0=[F(u),Io(u)];
% fun=@(x0)Alfa_HA_Lleida(x0,u);
% options = optimset('MaxFunEvals',100000,'MaxIter',100000);  % Optimització dels paràmetres de la funció 
% OUT(u,:)=fminsearch(fun,x0,options);
% end
