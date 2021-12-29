close all
clear all

%% Model SEIR per la GRIP %%

format long

v_anys = [2010,2011,2012,2013,2014,2015,2016,2017,2018,2019];

% Dades
[num] = xlsread('CasosR15_18.xlsx'); % llegeixo l'excel com a base de dades
HA = readtable('HA_BCN.xlsx','Sheet',2); % llegeixo l'excel com a taula


% Paràmetres
DeltaT = 1;             % days 
Tini = 0; Tfin = 364;   % days
TfinCR = 371;           % days

Temps = Tini:DeltaT:Tfin; Npassos = length(Temps);
TempsCR = Tini:7:TfinCR;

alfa_1 = 1.99907566206338e-06;
alfa_2 = 2.97494957788777e-06;
alfa_3 = 3.97828357556244e-06;

llindar_1 = 483;         % ind
llindar_2 = 1642;        % ind  
beta = 1/4;              % 1/dia
gamma = 1/7;             % 1/dia
N = 5e6;                 % ind    

F = [0.01880;0.01666;0.01697;0.01717;0.01697;0.01684;0.01745;0.01790;0.01713;0.01710]; % factors de reducció entre 2010 i 2018

% Infectats Inicials
Io = [9;68;50;43;50;47;44;24;48;42];

% Inicialitzem matrius per a guardar els infectats del model i de casos reals
save_I_mod = []; save_I_real = [];

% Colors plots
color = [0 0.4470 0.7410]; colorCR = [0.6350 0.0780 0.1840]; colorHA = [0.8500 0.3250 0.0980];

for i = 1:length(v_anys)    

    Any = v_anys(i);
    S = zeros(Npassos,1);             
    E = zeros(Npassos,1);
    I = zeros(Npassos,1);
    R = zeros(Npassos,1);
    Y = Any-2009;
    f = F(Y);
    
    % Casos Reals
    CasosReals = num(:,Y+1);
    % Humitat Absoluta
%     hum_abs = HA.(['HA',num2str(Any)]); % llegeixo els valors d'HA i faig moving average
%     ha_av = movmean(hum_abs,7);
ha_av=HA.(['HAM7',num2str(Any)]);

    % Valors Inicials
    S(1) = f*N;
    E(1) = 0;
    I(1) = Io(Y);
    R(1) = 0;
    
    for t=2:Npassos
        if I(t-1) < llindar_1
            S(t) = S(t-1)-(alfa_1*S(t-1)*I(t-1))*DeltaT;
            E(t) = E(t-1)+(alfa_1*S(t-1)*I(t-1)-beta*E(t-1))*DeltaT;
            I(t) = I(t-1)+(beta*E(t-1)-gamma*I(t-1))*DeltaT;
            R(t) = R(t-1)+(gamma*I(t-1))*DeltaT;
        end
        if I(t-1) >= llindar_1 && I(t-1) < llindar_2
            S(t) = S(t-1)-(alfa_2*S(t-1)*I(t-1))*DeltaT;
            E(t) = E(t-1)+(alfa_2*S(t-1)*I(t-1)-beta*E(t-1))*DeltaT;
            I(t) = I(t-1)+(beta*E(t-1)-gamma*I(t-1))*DeltaT;
            R(t) = R(t-1)+(gamma*I(t-1))*DeltaT;
        end
        if I(t-1) >= llindar_2
            S(t) = S(t-1)-(alfa_3*S(t-1)*I(t-1))*DeltaT;
            E(t) = E(t-1)+(alfa_3*S(t-1)*I(t-1)-beta*E(t-1))*DeltaT;
            I(t) = I(t-1)+(beta*E(t-1)-gamma*I(t-1))*DeltaT;
            R(t) = R(t-1)+(gamma*I(t-1))*DeltaT;  
        end
    end
    save_I_mod = [save_I_mod, I];                     
    save_I_real = [save_I_real, CasosReals];          
end


%% SEARCH FOR:

% !!! The incubation period of influenza is usually two days but can range from one to four days.    !!!
% !!! A person may pass virus from 1 day before symptoms start through 5–7 days after illness onset. !!!

% Inici i final del pic d'infectats
    % Determining the Period of Exposure
    % 1. Identify the peak of the outbreak, which is the time period then the largest number of cases occurred.
    [max_infected_mod, ind_M_mod] = max(save_I_mod); [max_infected_real, ind_M_real] = max(save_I_real);
    % 2. Count back from the peak, the average incubation period for disease. Note that date.
    av_incub_real = []; av_incub_mod = []; ind_av_mod = []; ind_av_real = [];
    for j = 1:length(v_anys)
        av_incub_mod(j) = save_I_mod(ind_M_mod(j) - 2,j); av_incub_real(j) = save_I_real(ind_M_real(j),j);
        ind_av_mod(j) = ind_M_mod(j) - 2; ind_av_real(j) = ind_M_real(j);
    end
    % 3. Identify the earliest case in the outbreak and count back the minimum incubation period. Note that date. 70 infectats pillo
    loc_mod = []; loc_real = []; earliest_case_mod = []; earliest_case_real = []; 
    for k = 1:length(v_anys)
        loc_mod(k) = find(save_I_mod(:,k)>70,1); loc_real(k) = find(save_I_real(:,k)>70,1); 
        earliest_case_mod(k) = save_I_mod(loc_mod(k) - 1,k);
        earliest_case_real(k) = save_I_real(loc_real(k),k);
    end
    % 4. Identify the last case in the outbreak and count back the maximum incubation period. Note that date.
    loc_mod_last = []; loc_real_last = []; latest_case_mod = []; latest_case_real = []; 
    for k = 1:length(v_anys)
        loc_mod_last(k) = find(save_I_mod(:,k)>70,1,'last'); loc_real_last(k) = find(save_I_real(:,k)>70,1,'last'); 
        latest_case_mod(k) = save_I_mod(loc_mod_last(k) - 4,k);
        latest_case_real(k) = save_I_real(loc_real_last(k) - 1,k);
    end
    % The range of dates identified in Step 2-4 represent the most likely period of exposure and ideally should 
    % fall within a few days of each other. If the three dates are widely separated, this indicates that the incubation 
    % period has a wide range or the outbreak is not a point source outbreak.

f = figure(1);
f.Position = [0 0 2400 1800];
for i = 1:length(v_anys)
    subplot(2,5,i)
    plot(Temps, save_I_mod(:,i), 'LineWidth', 1.2)
    hold on
    plot(TempsCR, save_I_real(:,i),'.r', 'MarkerSize',8)
    xline(TempsCR(loc_real_last(i)),'r--'); xline(loc_mod_last(i),'b--');
    xline(TempsCR(loc_real(i)),'r:','LineWidth',1.1); xline(loc_mod(i),'b:','LineWidth',1.15);
    xline(TempsCR(ind_av_real(i)),'r-.'); xline(ind_av_mod(i),'b-.');
    x = [TempsCR(loc_real_last(i)) TempsCR(loc_real(i)) TempsCR(loc_real(i)) TempsCR(loc_real_last(i))];
    y = [0 0 max(max_infected_real(i), max_infected_mod(i)) max(max_infected_real(i), max_infected_mod(i))];
    p = patch(x,y,'red'); alpha(p,.1);
    x2 = [loc_mod_last(i) loc_mod(i) loc_mod(i) loc_mod_last(i)];
    y2 = [0 0 max(max_infected_real(i), max_infected_mod(i)) max(max_infected_real(i), max_infected_mod(i))];
    p2 = patch(x2,y2,'blue'); alpha(p2,.05);
    hold off
    grid on
    legend('Model', 'Casos', Location='northwest')
    xlabel('Time (days)')
    ylabel('Infected individuals')
    title(['Period of exposure for ', num2str(v_anys(i))])
    ylim([0 15e3])
end

f2 = figure(2);
f2.Position = [0 0 2400 1800];
for i = 1:length(v_anys)
    subplot(2,5,i)
    plot(Temps(1:end-1), diff(save_I_mod(:,i))./diff(Temps'), 'LineWidth', 1.2)
    hold on
    plot(TempsCR(1:end-1), diff(save_I_real(:,i))./diff(TempsCR'),'r', 'LineWidth',1)
    hold off
    grid on
    legend('Model','Cases',Location='northwest')
    xlabel('Time (days)')
    ylabel('dI/dt')
    title(['Tax of growth ', num2str(v_anys(i))])
end

der_infect_mod = []; der_infect_real = [];

for l=1:length(v_anys)
    der_infect_mod(:,l) = diff(save_I_mod(:,l))./diff(Temps');
    der_infect_real(:,l) = diff(save_I_real(:,l))./diff(TempsCR');
end

% [pks,locs,w,p] = findpeaks(der_infect_real(:,1));
% 
% figure(3)
% plot(TempsCR(1:end-1), der_infect_real(:,1))
% hold on
% plot(TempsCR(locs), pks, '*m')
% hold off
% xlabel('Time (weeks)')
% ylabel('dI/dt with marked peaks')
% 
% figure(4)
% bar(p) % https://es.mathworks.com/help/signal/ref/findpeaks.html#buff2uu
% xlabel('Number of peak')
% ylabel('Prominence of the peaks')
% 
% figure(5)
% bar(w)
% xlabel('Number of peak')
% ylabel('Width of the peaks')

% vull marcar inici com primer pic amb prominència superior a mitjana abans del pic màxim, és a dir, quan la tendència comença una pujada important.
%%
[max_dinfdt_real, ind_M_didt_real] = max(der_infect_real);

loc_inici = []; infectats_inici = []; loc_inici_days = [];
for mm = 1:length(v_anys)
    [pks,locs,w,p] = findpeaks(der_infect_real(:,mm));
    when_max = find(locs == ind_M_didt_real(mm));
    threshold = mean(p(1:when_max-1));
    for extra = 1:when_max-1
        if abs(pks(extra) - pks(when_max)) < 50
            threshold = mean(p(1:extra-1));
        end
    end
    pic_inici =  find(p > threshold, 1);
    loc_inici(mm) = locs(pic_inici);
    infectats_inici(mm) = save_I_real(loc_inici(mm),mm);
    if infectats_inici(mm) < 50
        pic_inici =  find(p > threshold, 2);
        loc_inici(mm) = locs(pic_inici(2));
        infectats_inici(mm) = save_I_real(loc_inici(mm),mm);
    end
    loc_inici_days(mm) = loc_inici(mm)*7; % shiiiit
end

figure(3)
for i=1:length(v_anys)
    subplot(2,5,i)
    plot(TempsCR(loc_inici(i)),infectats_inici(i), 'bo')
    hold on
    plot(TempsCR, save_I_real(:,i),'.r', 'MarkerSize',8)
    plot(TempsCR(1:end-1), diff(save_I_real(:,i))./diff(TempsCR'),':k', 'MarkerSize',6)
    hold off
    title(['Year ', num2str(v_anys(i))])
    legend('Inicial peak', 'Cases', 'Tax of growth', Location='northwest')
end

%% Quin és el valor d'HA en la setmana en que l'epidèmia "comença" (comença a partir d'un nombre d'infectats a triar amb criteri propi).

humitats = []; humitats_7 = []; humitats_14 = []; humitats_2 = []; humitats_4 = [];
M_humitats = []; M_humitats_7 = []; M_humitats_14 = []; M_humitats_2 = []; M_humitats_4 = [];
setmanal_ha_av = [];

for i = 1:length(v_anys)    
    Any = v_anys(i);
    % Humitat Absoluta
    hum_av = HA.(['HAM7',num2str(Any)]); % llegeixo els valors d'HA i faig moving average
%     ha_av = movmean(hum_abs,7);

    ha_av_setm = []; ha_av_setm(1) = ha_av(1);
    for ss = 8:7:length(ha_av)
        ha_av_setm(ss) = mean(ha_av(ss-7:ss));
    end
    ha_av_setm(ha_av_setm==0) = [];
    ha_av_setm(54) = ha_av(end);

    humitats(i) = ha_av(loc_inici_days(i)); 
    humitats_7(i) = ha_av(loc_inici_days(i)-7);
    humitats_14(i) = ha_av(loc_inici_days(i)-14);
    humitats_2(i) = ha_av(loc_inici_days(i)-2);
    humitats_4(i) = ha_av(loc_inici_days(i)-4);
    
    M_humitats(i) = ha_av(ind_M_real(i)); 
    M_humitats_7(i) = ha_av(ind_M_real(i)-7);
    M_humitats_14(i) = ha_av(ind_M_real(i)-14);
    M_humitats_2(i) = ha_av(ind_M_real(i)-2);
    M_humitats_4(i) = ha_av(ind_M_real(i)-4);

    setmanal_ha_av(:,i) = ha_av_setm;
end

    figure
    subplot(2,3,1)
    scatter(humitats, infectats_inici)
    xlabel('HA initial day')
    ylabel('Infected people initial day')
    subplot(2,3,2)
    scatter(humitats_2, infectats_inici)
    xlabel('HA 2 days before')
    ylabel('Infected people initial day')
    subplot(2,3,3)
    scatter(humitats_4, infectats_inici)
    xlabel('HA 4 days before')
    ylabel('Infected people initial day')
    subplot(2,3,4)
    scatter(humitats_7, infectats_inici)
    xlabel('HA a week before')
    ylabel('Infected people initial day')
    subplot(2,3,6)
    scatter(humitats_14, infectats_inici)
    mdl = fitlm(humitats_14,infectats_inici);
    hold on
    plot(mdl)
    hold off
    xlabel('HA 2 weeks before')
    ylabel('Infected people initial day')




    figure
    subplot(2,3,1)
    scatter(M_humitats, max_infected_real)
    xlabel('HA initial day')
    ylabel('Max infected people')
    subplot(2,3,2)
    scatter(M_humitats_2, max_infected_real)
    xlabel('HA 2 days before')
    ylabel('Max infected people')
    subplot(2,3,3)
    scatter(M_humitats_4, max_infected_real)
    xlabel('HA 4 days before')
    ylabel('Max infected people')
    subplot(2,3,4)
    scatter(M_humitats_7, max_infected_real)
    xlabel('HA a week before')
    ylabel('Max infected people')
    subplot(2,3,6)
    scatter(M_humitats_14, max_infected_real)
    xlabel('HA 2 weeks before')
    ylabel('Max infected people')
%%
figure
for i = 1:length(v_anys)
    subplot(2,5,i)
    scatter(setmanal_ha_av(:,i),save_I_real(:,i))
    xlabel('HA (g/m^3)')
    ylabel('Infected people')
    title(['Year ', num2str(v_anys(i))])
    xlim([0 11])
end
