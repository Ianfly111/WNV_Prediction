clear; clc; close all

Case_NE_new = [1,1,6,5,11,16,15,17,20,15,18,13,7,1,0,1,1]; % start july 8 weekly, epiweek 27
Case_NE_cum = cumsum(Case_NE_new);
Pop_NE = 1000000;

Case_AR_new = [1,1,1,2,6,8,9,10,11,5,8,8,1,1,1,1,0,1,3,2,0,0,1,0,0,1]; % start epiweek 19
Case_AR_cum = cumsum(Case_AR_new);
Pop_AR = 7276000;

Case_CO_new = [1,1,1,2,10,25,38,68,92,86,83,79,55,42,19,13,6,4,0,3,0,0,0,0,0,0]; % start epiweek 24
Case_CO_cum = cumsum(Case_CO_new);
Case_CO_cum = Case_CO_cum(4:end);% start epiweek 27
Pop_CO = 3519582;
% Case_CO_cdc = [5,31,31,101,102,243,368,438,542,542,584];
% Case_CO = Case_CO_cdc;
% for i = 2:(length(Case_CO)-1)
%     Case_CO(i) = 0.5*Case_CO_cdc(i) + 0.25*(Case_CO_cdc(i-1) + Case_CO_cdc(i+1));
% end

Case = Case_CO_cum;
Pop = Pop_CO;
Nt = length(Case)-3;


y = Case(1:Nt)/Pop;
day = datetime('01-Jan-2023')+26*7;
Nens = 200;

d = 12;
C = [0 0 0 0 0 0 0 0 0 0 1 0];%; 0 0 0 0 1 1 1 0 0 0 0
q = [1e-4 1e-4 1e-4 1e-5 1e-3 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-2]/10;
x0 = zeros(d,Nens);
for i=1:Nens
    Eh = normrnd(0.00003,0.00001);
    Rh = 0;%60790*3/Pop
    Ch = 0;
    Sh = 1-Eh-Rh-Ch;
    
    Sm = normrnd(0.01,0.001);
    Em = normrnd(0.0003,0.0001);
    Im = normrnd(0.0001,0.00001); %normrnd(0,0.0000001);
    
    Sb = normrnd(0.01,0.002);
    Eb = normrnd(0.0003,0.0001);
    Ib = normrnd(0.0001,0.00001);
    Rb = normrnd(0.00001,0.000001);
    
    beta = 0.5;% CO 0.15
    
    x0(:,i) = [Sb, Eb, Ib, Rb, Sm, Em, Im, Sh, Eh, Rh, Ch, beta];%,BetaM
end
% **********************************************%

[No,Nt] = size(y);
X_pri = zeros(d,Nens,Nt);
M_pri = zeros(d,Nt);
X = zeros(d,Nens,Nt+1);
M = zeros(d,Nt);
P_pri = zeros(d,d,Nt);
P = zeros(d,d,Nt);
K = zeros(d,No,Nt);
X(:,:,1) = x0;

m = 3;
Fcst = cell(1,m);  
for k = 1:m
    Fcst{k} = zeros(d,Nt);
end

Nit = 500; 
Est =zeros(d,Nt);
Prt =cell(1,m);
for k = 1:m
    Prt{k} = zeros(Nit,d,Nt);
end
Prt_d = zeros(Nit,m);

for nn = 1:Nit
    for ii = 1:Nt
        V = diag(y(:,ii))/(3^ii);% with noise 
        Q = diag(q)/(ii);
   %prediction
       for j = 1:Nens
           X_pri(:,j,ii) = Forcast(X(:,j,ii));
       end 

       X_pri(:,:,ii) = max(0, X_pri(:,:,ii) + Q*randn(size(X_pri(:,:,ii))));
       M_pri(:,ii) = mean(X_pri(:,:,ii),2); % mean of priors
       LP_pri = 1/sqrt(Nens-1)*(X_pri(:,:,ii) - M_pri(:,ii)); % LP_pri*LP_pri' = P_pri
       P_pri(:,:,ii) = LP_pri*LP_pri';

       %ensemble observation noise
       E = V*randn(No,Nens);
       E = E - mean(E,2);
       Ve = 1/(Nens-1)*(E*E');
       Y = y(:,ii)+max(0,E);

       %filtering
       K(:,:,ii) = LP_pri*LP_pri'*C'*inv(C*(LP_pri*LP_pri')*C'+Ve); %Gain
       for j = 1:Nens
           X(:,j,ii+1) = max(0,X_pri(:,j,ii) + K(:,:,ii)*(Y(:,j)-C*X_pri(:,j,ii)));
       end
       M(:,ii) = 0.5*max(0,M_pri(:,ii) + K(:,:,ii)*(y(:,ii)-C*M_pri(:,ii)))+0.5*max(0,mean(X(:,:,ii+1),2));
       LP = 1/sqrt(Nens-1)*(X(:,:,ii+1) - M(:,ii));
       P(:,:,ii) = LP*LP';

       Fcst{1}(:,ii) = max(0,Forcast(M(:,ii)));% +K(:,:,ii)*(y(:,ii)-C*M_pri(:,ii))
       Fcst{1}(11,ii) = max(Fcst{1}(11,ii),M(11,ii));
       if m>1
           for k = 2:m
               Fcst{k}(:,ii) = max(0,Forcast(Fcst{k-1}(:,ii)));%+K(:,:,ii)*(y(:,ii)-C*M_pri(:,ii))/k
               Fcst{k}(11,ii) = max(Fcst{k}(11,ii),Fcst{k-1}(11,ii));
           end
       end 
       
    end
    
    Est = Est+M/Nit;
    for k = 1:m
        Prt{k}(nn,:,:) = Fcst{k};
    end
    
    for k = 1:m
        Prt_d(nn,k) = Fcst{k}(11,ii);
    end
end

Days = day + (1:Nt)*7;
figure
subplot(4,1,1)
plot(Days,Est(1,:));
legend('estimated Sb')
xlabel 'Time';
ylabel 'Sb'

subplot(4,1,2)
plot(Days,Est(2,:));
legend('estimated Eb')
xlabel 'Time';
ylabel 'Eb'

subplot(4,1,3)
plot(Days,Est(3,:));
legend('estimated Ib')
xlabel 'Time';
ylabel 'Ib'

subplot(4,1,4)
plot(Days,Est(4,:));
legend('estimated Rb')
xlabel 'Time';
ylabel 'Rb'

figure
subplot(3,1,1)
plot(Days,Est(5,:));
legend('estimated Sm')
xlabel 'Time';
ylabel 'Sm'

subplot(3,1,2)
plot(Days,Est(6,:));
legend('estimated Em')
xlabel 'Time';
ylabel 'Em'

subplot(3,1,3)
plot(Days,Est(7,:));
legend('estimated Im')
xlabel 'Time';
ylabel 'Im'


figure
subplot(4,1,1)
plot(Days,Est(8,:));
legend('estimated S_{h}')
xlabel 'Time';
ylabel 'Sh'

subplot(4,1,2)
plot(Days,Est(9,:));
legend('estimated E_{h}')
xlabel 'Time';
ylabel 'Eh'

subplot(4,1,3)
plot(Days,Est(10,:));
legend('estimated R_{h}')
xlabel 'Time';
ylabel 'Rh'

subplot(4,1,4)
plot(Days,Est(11,:));
legend('estimated C_{h}')
xlabel 'Time';
ylabel 'Ch'

figure
plot(Days,Est(12,:));
legend('estimated beta')
xlabel 'Time';
ylabel 'infection rate'

figure
boxplot(round(Prt_d*Pop),'Symbol',' ');
xlabel 'Time (weeks)';
ylabel 'Distribution'

p1 = round(mean(Prt{1}(:,11,:)*Pop));
p2 = round(mean(Prt{2}(:,11,:)*Pop));
p3 = round(mean(Prt{3}(:,11,:)*Pop));
figure
plot(1:(Nt+2),Case(2:Nt+3), 1:Nt, reshape(p1,[1,Nt]), 2:(Nt+1), reshape(p2,[1,Nt]), 3:(Nt+2), reshape(p3,[1,Nt]))
legend('Reported cases', 'Ave 1-week', 'Ave 2-week', 'Ave 3-week')
xlabel 'Time (weeks)';
ylabel 'Symptomatic cases'

Daysbox = zeros(1,Nt*Nit);
pbox1 = zeros(1,Nt*Nit);
pbox2 = zeros(1,Nt*Nit);
pbox3 = zeros(1,Nt*Nit);
pnbox1 = zeros(1,Nt*Nit);
pnbox2 = zeros(1,Nt*Nit);
pnbox3 = zeros(1,Nt*Nit);
for i=1:Nt
    t1 = Nit*(i-1)+1;
    t2 = Nit*i;
    Daysbox(t1:t2) = i*7;
    pbox1(t1:t2) = round(Prt{1}(:,11,i)*Pop);
    pbox2(t1:t2) = round(Prt{2}(:,11,i)*Pop);
    pbox3(t1:t2) = round(Prt{3}(:,11,i)*Pop);
    pnbox1(t1:t2) = round(max(0,Prt{1}(:,11,i)-y(i))*Pop);
    pnbox2(t1:t2) = round(max(0,Prt{2}(:,11,i)-Prt{1}(:,11,i))*Pop);
    pnbox3(t1:t2) = round(max(0,Prt{3}(:,11,i)-Prt{2}(:,11,i))*Pop);
end

figure
plot(Case(2:Nt+1),'*')
hold on
boxplot(pbox1,day+7+Daysbox, 'symbol','')
hold off
xlabel 'Time';
ylabel 'One-week predicted symptomatic cases'

figure
plot(Case(3:Nt+2),'*')
hold on
boxplot(pbox2, day+14+Daysbox, 'symbol','')
hold off
xlabel 'Time';
ylabel 'Two-week predicted symptomatic cases'

figure
plot(Case(4:Nt+3),'*')
hold on
boxplot(pbox3, day+21+Daysbox, 'symbol','')
hold off
xlabel 'Time';
ylabel 'Three-week predicted symptomatic cases'

Case_n = diff(Case);
figure
plot(Case_n(1:Nt),'*')
hold on
boxplot(pnbox1,day+7+Daysbox, 'symbol','')
hold off
xlabel 'Time';
ylabel 'One-week predicted new symptomatic cases'

figure
plot(Case_n(2:Nt+1),'*')
hold on
boxplot(pnbox2, day+14+Daysbox, 'symbol','')
hold off
xlabel 'Time';
ylabel 'Two-week predicted new symptomatic cases'

figure
plot(Case_n(3:Nt+2),'*')
hold on
boxplot(pnbox3, day+21+Daysbox, 'symbol','')
hold off
xlabel 'Time';
ylabel 'Three-week predicted new symptomatic cases'

bin = zeros(2,Nt);
for i = 1:Nt+2
    if Case_n(i) <=5
        bl = 0;
        bh = 5;
    elseif Case_n(i) <=10
        bl = 5;
        bh = 10;
    elseif Case_n(i) <=15
        bl = 10;
        bh = 15;
    elseif Case_n(i) <=20
        bl = 15;
        bh = 20;
    elseif Case_n(i) <=25
        bl = 20;
        bh = 25;
    elseif Case_n(i) <=30
        bl = 25;
        bh = 30;
    elseif Case_n(i) <=35
        bl = 30;
        bh = 40;
    elseif Case_n(i) <=45
        bl = 40;
        bh = 45;
    elseif Case_n(i) <=50
        bl = 45;
        bh = 50;
    elseif Case_n(i) <=100
        bl = 50;
        bh = 100;
    else
        bl = 100;
        bh = inf;    
    end 
    bin(:,i) = [bl;bh];
end
lsn1 = zeros(1,Nt);
lsn2 = zeros(1,Nt);
lsn3 = zeros(1,Nt);
rmse1 = zeros(1,Nt);
rmse2 = zeros(1,Nt);
rmse3 = zeros(1,Nt);
m = -10;
for i = 1:Nt
    t1 = Nit*(i-1)+1;
    t2 = Nit*i;
    rmse1(i) = sqrt(mean((pbox1(t1:t2)-Case(i+1)).^2));
    rmse2(i) = sqrt(mean((pbox2(t1:t2)-Case(i+2)).^2));
    rmse3(i) = sqrt(mean((pbox3(t1:t2)-Case(i+3)).^2));
    lsn1(i) = max(m,log(sum((pnbox1(t1:t2)>bin(1,i)).*(pnbox1(t1:t2)<=bin(2,i)))/Nit));
    lsn2(i) = max(m,log(sum((pnbox2(t1:t2)>bin(1,i+1)).*(pnbox2(t1:t2)<=bin(2,i+1)))/Nit));
    lsn3(i) = max(m,log(sum((pnbox3(t1:t2)>bin(1,i+2)).*(pnbox3(t1:t2)<=bin(2,i+2)))/Nit));
end

figure
plot(1:Nt,lsn1,1:Nt,lsn2,1:Nt,lsn3)
legend('One-week','Two-week','Three-week')
ylabel('Logarithmic score')
ylim([0.5*m,0.5])
xlabel('Week')

figure
plot(1:Nt,rmse1,1:Nt,rmse2,1:Nt,rmse3)
legend('One-week','Two-week','Three-week')
ylabel('Rooted mean square error')
%ylim([0.5*m,0.5])
xlabel('Week')


% RTPreR = zeros(m,Nt);
% RTPreI = zeros(m,Nt);
% Error = zeros(m,Nt);
% for R = 1:m
%    RTPreR(R,:) = Prt{R}(4,:); 
%    RTPreI(R,:) = Prt{R}(3,:);
%    Error(R,:) = abs(RTPreR(R,:)-CumCase(1+R:30+R))./CumCase(1+R:30+R);
% end
% t = 1:Nt;
% figure
% subplot(2,1,1)
% plot(1:length(CumCase),CumCase, t+1,RTPreR(1,:), t+2,RTPreR(2,:), t+3,RTPreR(3,:), t+4,RTPreR(4,:), t+5,RTPreR(5,:));
% legend('R_{h}','R1_{h}','R2_{h}','R3_{h}','R4_{h}','R5_{h}')
% xlabel 'Time (weeks)';
% ylabel 'Number of incidence'
% 
% subplot(2,1,2)
% plot(t+1,RTPreI(1,:), t+2,RTPreI(2,:), t+3,RTPreI(3,:), t+4,RTPreI(4,:), t+5,RTPreR(5,:));
% ylim([0, 0.001])
% legend('I1_{h}','I2_{h}','I3_{h}','I4_{h}','I5_{h}')
% xlabel 'Time (weeks)';
% ylabel 'Predicted infectious people'
% 
% figure
% plot(t+1,Error(1,:), t+2,Error(2,:), t+3,Error(3,:), t+4,Error(4,:), t+5,Error(5,:));
% ylim([0,2])
% legend('E1_{h}','E2_{h}','E3_{h}','E4_{h}','E5_{h}')
% xlabel 'Time (weeks)';
% ylabel 'Relative error'


function P=Forcast(SEIRSEISEIR)
% sys state-space function
    Sb = SEIRSEISEIR(1); Eb = SEIRSEISEIR(2); Ib = SEIRSEISEIR(3); Rb = SEIRSEISEIR(4);
    Sm = SEIRSEISEIR(5); Em = SEIRSEISEIR(6); Im = SEIRSEISEIR(7);
    Sh = SEIRSEISEIR(8); Eh = SEIRSEISEIR(9); Rh = SEIRSEISEIR(10); Ch = SEIRSEISEIR(11);
    
    beta = max(0, SEIRSEISEIR(12));
    
    
    vb = 0.05*normrnd(1,0.2); mub = 0.05*normrnd(1,0.3); betamb = 5*beta*normrnd(1,0.3); deltab = 2/3*normrnd(1,0.2); gammab = 1/3*normrnd(1,0.2); 
    vm = 0.3*normrnd(1,0.2); mum = 1/6*normrnd(1,0.2); betabm = 10*beta*normrnd(1,0.2); deltam = 1/2*normrnd(1,0.2);
    vh = 0; muh = 0; betamh = beta; deltah = 0.6*normrnd(1,0.2); gammah = 0.06*normrnd(1,0.2);
    
    Nb = Sb + Eb + Ib + Rb;
    Nm = Sm + Em + Im;
    
    Sb1= vb*Nb - (mub + betamb*Im)*Sb + Sb; 
    Eb1= betamb*Im*Sb - (mub+deltab)*Eb + Eb;
    Ib1= deltab*Eb - (mub+gammab)*Ib + Ib; 
    Rb1= gammab*Ib - mub*Rb + Rb;
    
    Sm1= vm*Nm -(mum + betabm*Ib)*Sm +Sm;  
    Em1= betabm*Ib*Sm - (mum + deltam)*Em + Em;
    Im1= deltam*Em - mum*Im +Im;
    
    Sh1= vh - (muh + betamh*Im)*Sh + Sh;
    Eh1= betamh*Im*Sh - (deltah + gammah + muh)*Eh + Eh;
    Rh1= deltah*Eh - muh*Rh +Rh;
    Ch1= gammah*Eh - muh*Ch +Ch;
    
    P = [Sb1; Eb1; Ib1; Rb1; Sm1; Em1; Im1; Sh1; Eh1; Rh1; Ch1; beta*normrnd(1,0.2)];
end