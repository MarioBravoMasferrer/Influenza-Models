%% Gràfics dades reals
close all
[num]=xlsread('CasosR15_18.xlsx');

for v=1:10
d1=cell(52,1);
for p=23:52
    d1(p-22)={num2str(p)};
end
for p=1:22
    d1(p+30)={num2str(p)};
end
f1=figure(1);
f1.Position = [100 10 1100 1200];
DarkBlue=[0 0.2 0.5];
yellow = 1/255*[240,180,15];
l=v;
if v>8
    l=v+1;
end
subplot(3,4,l)
plot(num(1:52,v+1),'.','MarkerSize',7,'Color',yellow)
set(gca,'xtick',1:6:52);
set(gca,'xticklabel',d1(1:6:end))
hold on
title([num2str(2009+v),'-',num2str(2009+v+1)])
xlabel("Setmanes de l'any")
if v==1 || v==5 || v==9
ylabel('Individus Infectats')
end
xlim([0 54])
ylim([0 13000])
set(gca,'ytick',0:2000:13000);
hold off
end

d1=cell(52,1);
for p=23:52
    d1(p-22)={num2str(p)};
end
for p=1:22
    d1(p+30)={num2str(p)};
end
pink=[0.7 0.2 0.8];
f2=figure(20);
plot(num(1:52,2),'.-','MarkerSize',10,'Color',DarkBlue)
hold on
plot(num(1:52,3),'.-','MarkerSize',10)
plot(num(1:52,4),'.-','MarkerSize',10)
plot(num(1:52,5),'.-','MarkerSize',10)
plot(num(1:52,6),'.-','MarkerSize',10)
plot(num(1:52,7),'.-','MarkerSize',10)
plot(num(1:52,8),'.-','MarkerSize',10)
plot(num(1:52,9),'.-','MarkerSize',10)
plot(num(1:52,10),'.-','MarkerSize',10,'Color','r')
plot(num(1:52,11),'.-','MarkerSize',10,'Color',pink)
set(gca,'xtick',1:1:52);
set(gca,'xticklabel',d1(1:1:end))
hold on
a=cell(10,1);
for v=1:10
   a(v)={[num2str(2009+v),'-',num2str(2009+v+1)]}; 
end
legend(a)
xlabel("Setmanes de l'any")
ylabel('Individus Infectats')
xlim([17 52])
ylim([0 13000])
set(gca,'ytick',0:1000:13000);
yyaxis right
yticks([])

hold off

for v=2:10
o=1;
for i=1:52
    if num(i,v+1)>800 && o==1
      loc_900=i;
      o=0;
    end
end
for i=1:52
    if num(i,v+1)>6700 && o==0
      loc_6400=i;
      o=1;
    end
end
f3=figure(3);
f3.Position = [100 10 1100 1200];
DarkBlue=[0 0.2 0.5];
yellow = 1/255*[240,180,15];
l=v;
subplot(3,3,l-1)
plot(num(1:52,v+1),'.','MarkerSize',7,'Color',yellow)
set(gca,'xtick',1:6:52);
set(gca,'xticklabel',d1(1:6:end))
orange=[1 0.5 0];
xline(loc_900,'Color',orange,'LineWidth',1.2)
xline(loc_6400,'Color','r','LineWidth',1.2)
text(loc_900,10000,'\leftarrow')
text(loc_6400-3,10000,'\rightarrow')
text(loc_900-12,10000,[num2str(loc_6400-loc_900),' Setmanes'],'FontSize',7)
hold on
title([num2str(2009+v),'-',num2str(2009+v+1)])
if v>7
xlabel("Setmanes de l'any")
end
if v==2 || v==5 || v==8
ylabel('Individus Infectats')
end
xlim([0 54])
ylim([0 13000])
set(gca,'ytick',0:2000:13000);
if v==2
legend('Casos Reals','800 infectats', '6700 infectats')
end
hold off
end

