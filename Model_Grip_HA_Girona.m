% function [error]=Alfa_HA_Girona(x0,u)
format long
close all
%% ALFA(HA)
%% Paràmetres
v_anys = [2011,2012,2013,2014,2015,2016,2017,2018,2019];
% v_anys=v_anys(u);
[num]=xlsread('Grip Girona.xlsx');                     
HA=xlsread('HA_Girona');        
HA=HA(1:end,1:9);

DeltaT = 1;           %dia
Tini=1;               %dies
Tfin=357;             %dies
beta = 1/4;           %1/dia
gamma = 1/7;          %1/dia
N=888467;            %ind   
NB=5e6;
alfas=[];
F=[0.0149236105238601;0.0169295839156001;0.0169888867322746;0.0179494560593738;0.0165644085798682;0.0177732047358208;0.0182873407065669;0.0184661827419294;0.0170005620292491];
Io=round([14.8109703769275;0.969778717058019;1.03997821914555;1.81849340304685;0.64735057304347;4.75481443808862;2.48612233803078;1.04498484464146;1.04608626319666]);
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
% % Gràfiques
Ini23=[7,6,4,3,2,1,6,5,4,3];
Blue = 1/255*[0,150,200];
v=i;
if i>8
    v=i+1;
end

subplot(3,3,i)
hold on
Infectats=plot(Temps,I(:),'k','LineWidth',1.2);
hold on
plot(TempsCR,CasosReals,'.','Color',Blue,'MarkerSize',8)
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
% OUT=zeros(9,2);
% for u=1:9
% F=[0.0157905682864913,0.0153164471970013,0.0179776884551847,0.0168353145218735,0.0170066415007447,0.0161399935283330,0.0183365014233892,0.0185843401086952,0.0193544254824303,0.0179297124222972];
% Io=[6,6,1,1,3,3,3,1,1,1];
% x0=[F(u),Io(u)];
% fun=@(x0)Alfa_HA_Girona(x0,u);
% options = optimset('MaxFunEvals',100000,'MaxIter',100000);  % Optimització dels paràmetres de la funció 
% OUT(u,:)=fminsearch(fun,x0,options);
% end
