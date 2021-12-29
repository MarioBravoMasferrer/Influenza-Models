%% Model SEIR per la GRIP %%
%% Parametritzada amb dades de la regió sanitària de Barcelona %%
%% EN AQUETS MODEL S'HAN PROVAT I BUSCAT VARIES REPSOSTES (Optimitzar les alfes i llindars, variació de f,...)

close all
% v_anys = linspace(2010,2019,10);
v_anys = [2010,2011,2012,2013,2014,2015,2016,2017,2018,2019];
[num]=xlsread('CasosR15_18.xlsx');
R2=[];R2G=[];
p_val=[];p_valG=[];

% Paràmetres
DeltaT = 1;           %dia
Tini=0;               %dies
Tfin=357;             %dies



alfa_1=2e-6;
alfa_2=3e-6;
alfa_3=4.3e-6;



llindar1=500;
llindar2=1500;

beta = 1/4;           %1/dia
gamma = 1/7;          %1/dia
N=5e6;            %ind    



% Vectors
Temps=Tini:DeltaT:Tfin;
Npassos=length(Temps);
F=[0.01722;0.01666;0.01697;0.01717;0.01697;0.01684;0.01745;0.01790;0.01717;0.01710];
TempsCR=0:7:364;
Error=[];
IM7=zeros(52,length(v_anys));
error=0;


% Infectats Inicials
Io = [55;68;50;43;50;47;44;24;48;42];



for i = 1:length(v_anys)    
Any = v_anys(i);
S=zeros(Npassos,1);
E=zeros(Npassos,1);
I=zeros(Npassos,1);
R=zeros(Npassos,1);
Y=Any-2009;
f=F(Y);


% Casos Reals
CasosReals=num(:,Y+1);
% Io=CasosReals(1);
% Valors Inicials
S(1)=f*N;
E(1)=0;
I(1)=Io(Y);
R(1)=0;

% MatrixI(1,i)=Io(Y);


for t=2:Npassos
    if I(t-1)<llindar1
    S(t)=S(t-1)-(alfa_1*S(t-1)*I(t-1))*DeltaT;
    E(t)=E(t-1)+(alfa_1*S(t-1)*I(t-1)-beta*E(t-1))*DeltaT;
    I(t)=I(t-1)+(beta*E(t-1)-gamma*I(t-1))*DeltaT;
    R(t)=R(t-1)+(gamma*I(t-1))*DeltaT;
    end
    if I(t-1)>=llindar1 && I(t-1)<llindar2
        S(t)=S(t-1)-(alfa_2*S(t-1)*I(t-1))*DeltaT;
        E(t)=E(t-1)+(alfa_2*S(t-1)*I(t-1)-beta*E(t-1))*DeltaT;
        I(t)=I(t-1)+(beta*E(t-1)-gamma*I(t-1))*DeltaT;
        R(t)=R(t-1)+(gamma*I(t-1))*DeltaT;
    end
    if I(t-1) >= llindar2
        S(t)=S(t-1)-(alfa_3*S(t-1)*I(t-1))*DeltaT;
        E(t)=E(t-1)+(alfa_3*S(t-1)*I(t-1)-beta*E(t-1))*DeltaT;
        I(t)=I(t-1)+(beta*E(t-1)-gamma*I(t-1))*DeltaT;
        R(t)=R(t-1)+(gamma*I(t-1))*DeltaT;  
    end

end
[r,p]=corrcoef(I(1:7:end),CasosReals(1:52));
R2(i)=r(1,2);
p_val(i)=p(1,2);


% 
%% Gràfiques 
% 
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
yellow = 1/255*[240,180,15];
v=i;
if i>8
    v=i+1;
end
% d1 = datetime(Any,6,Ini23(Y)) + caldays(1:length(I));
subplot(3,4,v)
hold on
Infectats=plot(Temps,I(:),'k','LineWidth',1.2);
hold on
plot(TempsCR,CasosReals,'.','Color',yellow,'MarkerSize',8)
set(gca,'xtick',1:28:length(I));
set(gca,'xticklabel',d1(1:4:end))
title([num2str(Any),'-',num2str(Any+1)])
xlabel("Setmanes de l'any");
if i==1 || i==5 || i==9
ylabel('Individus Infectats')
end
xlim([0 364])
ylim([0 14000])
set(gca,'ytick',0:2000:14000);
text(10,13000,['R2= ',num2str(R2(i))],'FontSize',7)
text(10,11500,['p-value= ',num2str(p_val(i))],'FontSize',7)
if i==1
legend('Model','Casos Reals','location','northwest');
end
hold off


% f2 = figure(2);
% f2.Position = [100 10 1100 1200];
% % Gràfiques
% subplot(3,4,i)
% hold on
% [~, idx] = max(I);
% [~, idx2] = max(CasosReals);
% plot(Temps(1:idx),I(1:idx));
% hold on
% plot(TempsCR(1:idx2),CasosReals(1:idx2),'.','MarkerSize',10)
% title(['Any ',num2str(Any),''])
% xlabel('Temps (dies)');
% ylabel('Infectats (ind)');
% ylim([0 1800])
% xlim([0 220])
% legend('Model','Casos Reals','location','northwest');
% hold off




%% Error



l=0;
Mid=floor(length(Temps(1:idx-l))/7);
new_I_matrix=zeros(Mid,length(v_anys));


len=floor(length(Temps(1:idx-l))/7);
new_I_matrix(1:len,i)=I(1:7:floor(length(Temps(1:idx-l)))-6);

for j=1:Mid
    error = error + (new_I_matrix(j,i)-CasosReals(j))^2;
end
Error(i)=error;

% [~, idx] = max(I);
[r,p]=corrcoef(I(1:7:floor(length(Temps(1:idx)))-6),CasosReals(1:idx/7));
R2G(i)=r(1,2);
p_valG(i)=p(1,2);

 end
error=sum(Error(:));



