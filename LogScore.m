clc; clear; close all
load('C:\Users\cyi\Desktop\Dengue CB\NewCaseCB.mat');
L = length(Case);
bin = zeros(2,L-41);
for i = 42:L
    if Case(i) == 0
        bl = 0;
        bh = 0;
    elseif Case(i) <=20
        bl = 0;
        bh = 20;
    elseif Case(i) <=40
        bl = 20;
        bh = 40;
    elseif Case(i) <=60
        bl = 40;
        bh = 60;
    elseif Case(i) <=80
        bl = 60;
        bh = 80;
    elseif Case(i) <=100
        bl = 80;
        bh = 100;
    elseif Case(i) <=150
        bl = 100;
        bh = 150;
    elseif Case(i) <=200
        bl = 150;
        bh = 200;
    elseif Case(i) <=250
        bl = 200;
        bh = 250;
    elseif Case(i) <=300
        bl = 250;
        bh = 300;
    elseif Case(i) <=350
        bl = 300;
        bh = 350;
    elseif Case(i) <=400
        bl = 350;
        bh = 400;
    elseif Case(i) <=450
        bl = 400;
        bh = 450;
    elseif Case(i) <=500
        bl = 450;
        bh = 500;
    else
        bl = 500;
        bh = inf;    
    end 
    bin(:,i-41) = [bl;bh];
end
ls_kf = zeros(1,length(bin));
ls_pf = zeros(1,length(bin));
ls_se = zeros(1,length(bin));
m = -5;
for i = 41:L-1
    week = mod(i,52);
    if week == 0
        week = 52;
    end
    kfm = strcat('kf_',num2str(week),'.mat');
    pfm = strcat('pf_',num2str(week),'.mat');
    sem = strcat('se_',num2str(week),'.mat');
    load(kfm);
    load(pfm);
    load(sem); 
    se(:,1) = round(0.4*kf(:,1)+0.6*pfn(:,1));
    save(sem,'se')
    ls_kf(i-40) = max(m,log(sum((kf(:,1)>bin(1,i-40)).*(kf(:,1)<=bin(2,i-40)))/500));
    ls_pf(i-40) = max(m,log(sum((pfn(:,1)>bin(1,i-40)).*(pfn(:,1)<=bin(2,i-40)))/500));
    ls_se(i-40) = max(m,log(sum((se(:,1)>bin(1,i-40)).*(se(:,1)<=bin(2,i-40)))/500));
end 

mean(ls_kf)
mean(ls_pf)
mean(ls_se)

ew = 42:L;
epiw = cell(size(ew));
for i = 1:length(ew)
    if ew(i) > 52
        w = ew(i) - 52;
        y = '23-';
    else
        w = ew(i);
        y = '22-';
    end
    epiw{i} = strcat(y,num2str(w));
end
epiw = categorical(epiw);

plot(epiw,ls_kf,'s',epiw,ls_pf,'*',epiw,ls_se,'o',epiw,zeros())
legend('Kalman Filter','Particle Filter','Super Ensemble')
ylabel('Logarithmic score')
ylim([1.3*m,0.5])
xlabel('Epiweek')