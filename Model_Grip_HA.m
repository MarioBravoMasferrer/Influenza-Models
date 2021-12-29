format long
close all

%% ALFA(HA)
%% Paràmetres
v_anys = [2010,2011,2012,2013,2014,2015,2016,2017,2018,2019];

[num]=xlsread('CasosR15_18.xlsx');                     
HA=xlsread('HA_BCN');        
HA=HA(1:end,:);

DeltaT = 1;           %dia
Tini=1;               %dies
Tfin=357;             %dies
beta = 1/4;           %1/dia
gamma = 1/7;          %1/dia
N=5e6;            %ind   

F=[0.0157905682864913,0.0153164471970013,0.0179776884551847,0.0168353145218735,0.0170066415007447,0.0161399935283330,0.0183365014233892,0.0185843401086952,0.0193544254824303,0.0179297124222972];  % Error 1.8e-8
Io=round([86.5677108089879,76.1353112729644,2.78364165584862,16.6338446284347,46.8231764506672,40.6135955050504,39.7313976521193,16.4180257271306,3.91231285302538,5.39172980254877]);
Temps=Tini:DeltaT:Tfin;
Npassos=length(Temps);
TempsCR=Tini:7:Tfin;
Error=zeros(1:length(v_anys));
error=0;
llindHA=[11.9478332519531,6.42500152587891];

Evol_CC=[10471354323.8632;4309046778.03316;5128191264.78385;3380671169.74687;2337501476.72572;1863965945.21516;1774374119.88371;873800985.754057;633628567.339877;576922496.482524;585389203.182084;417713198.432450;365675140.053636;375588744.682748;375588744.682748;375705294.582882;364698553.725737;364698553.725737];
HAvals=zeros(51,21);
R2G=[];p_valG=[];
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
HAvals(:,i*2)=HA(1:7:51*7);
HAvals(:,i*2+1)=CasosReals;

% Valors Inicials
S(1)=f*N;
E(1)=0;
I(1)=Io(Y);
R(1)=0;
alfa(1)=1.7e-06;


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


    S(t)=S(t-1)-(alfa(t-1)*S(t-1)*I(t-1))*DeltaT;
    E(t)=E(t-1)+(alfa(t-1)*S(t-1)*I(t-1)-beta*E(t-1))*DeltaT;
    I(t)=I(t-1)+(beta*E(t-1)-gamma*I(t-1))*DeltaT;
    R(t)=R(t-1)+(gamma*I(t-1))*DeltaT;   


end
[~, idx] = max(I);
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
ylim([0 16000])
set(gca,'ytick',0:2000:16000);
text(10,13500,['R2= ',num2str(R2(i))],'FontSize',7)
text(10,15000,['p-value= ',num2str(p_val(i))],'FontSize',7)
if i==1
legend('Model','Casos Reals','location','northwest');
end


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
error=sum(Error(:));




