clear; clc; close all
%load Case_GX;

Case_NE = [1,3,7,13,22,39,50,67,84,104,118,130,142];
Pop_NE = 1000000;

Case_CO_cdc = [5,31,31,101,102,243,368,438,542,542,584];
Case_CO = Case_CO_cdc;
for i = 2:(length(Case_CO)-1)
    Case_CO(i) = 0.5*Case_CO_cdc(i) + 0.25*(Case_CO_cdc(i-1) + Case_CO_cdc(i+1));
end

Case_CO_new = [1,0,0,1,0,1,0,2,0,10,24,38,68,89,86,83,78,55,42,19,13,6,4,0,3,0];
Case_CO_cum = cumsum(Case_CO_new);
Pop_CO = 3519582;


y = Case_CO_cum(9:20)/Pop_CO;% 9--week 27
day = datetime('01-Jan-2023')+26*7;
Nens = 200;

d = 12;
C = [0 0 0 0 0 0 0 0 0 0 1 0];%; 0 0 0 0 1 1 1 0 0 0 0
q = [1e-4 1e-4 1e-4 1e-5 1e-3 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-2]/10;
x0 = zeros(d,Nens);
for i=1:Nens
    Eh = normrnd(0.00003,0.00001);
    Rh = 60790*3/Pop_CO;%0
    Ch = 0;
    Sh = 1-Eh-Rh-Ch;
    
    Sm = normrnd(0.01,0.001);
    Em = normrnd(0.0003,0.0001);
    Im = normrnd(0.0001,0.00001); %normrnd(0,0.0000001);
    
    Sb = normrnd(0.01,0.002);
    Eb = normrnd(0.0003,0.0001);
    Ib = normrnd(0.0001,0.00001);
    Rb = normrnd(0.00001,0.000001);
    
    beta = 0.15;
    
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
       if m>1
           for k = 2:m
               Fcst{k}(:,ii) = max(0,Forcast(Fcst{k-1}(:,ii)));%+K(:,:,ii)*(y(:,ii)-C*M_pri(:,ii))/k
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
boxplot(round(Prt_d*Pop_CO),'Symbol',' ');
xlabel 'Time (weeks)';
ylabel 'Distribution'

p1 = round(mean(Prt{1}(:,11,:)*Pop_CO));
p2 = round(mean(Prt{2}(:,11,:)*Pop_CO));
p3 = round(mean(Prt{3}(:,11,:)*Pop_CO));
figure
plot(1:(Nt+2),Case_CO_cum(10:23), 1:Nt, reshape(p1,[1,Nt]), 2:(Nt+1), reshape(p2,[1,Nt]), 3:(Nt+2), reshape(p3,[1,Nt]))
legend('Reported cases', 'Ave 1-week', 'Ave 2-week', 'Ave 3-week')
xlabel 'Time (weeks)';
ylabel 'Symptomatic cases'

Daysbox = zeros(1,Nt*Nit);
pbox1 = zeros(1,Nt*Nit);
pbox2 = zeros(1,Nt*Nit);
pbox3 = zeros(1,Nt*Nit);
for i=1:Nt
    t1 = Nit*(i-1)+1;
    t2 = Nit*i;
    Daysbox(t1:t2) = i*7;
    pbox1(t1:t2) = round(Prt{1}(:,11,i)*Pop_CO);
    pbox2(t1:t2) = round(Prt{2}(:,11,i)*Pop_CO);
    pbox3(t1:t2) = round(Prt{3}(:,11,i)*Pop_CO);
end

figure
plot(Case_CO_cum(10:21),'*')
hold on
boxplot(pbox1,day+7+Daysbox, 'symbol','')
hold off
xlabel 'Time';
ylabel 'One-week predicted symptomatic cases'

figure
plot(Case_CO_cum(11:22),'*')
hold on
boxplot(pbox2, day+14+Daysbox, 'symbol','')
hold off
xlabel 'Time';
ylabel 'Two-week predicted symptomatic cases'

figure
plot(Case_CO_cum(12:23),'*')
hold on
boxplot(pbox3, day+21+Daysbox, 'symbol','')
hold off
xlabel 'Time';
ylabel 'Three-week predicted symptomatic cases'
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

% 
% function P=Forcast(SEIRSEISEIR)
% % sys state-space function
%     Sb = SEIRSEISEIR(1); Eb = SEIRSEISEIR(2); Ib = SEIRSEISEIR(3); Rb = SEIRSEISEIR(4);
%     Sm = SEIRSEISEIR(5); Em = SEIRSEISEIR(6); Im = SEIRSEISEIR(7);
%     Sh = SEIRSEISEIR(8); Eh = SEIRSEISEIR(9); Rh = SEIRSEISEIR(10); Ch = SEIRSEISEIR(11);
%     
%     beta = max(0, SEIRSEISEIR(12));
%     
%     
%     vb = 0.05*normrnd(1,0.2); mub = 0.05; betamb = 5*beta; deltab = 2/3; gammab = 1/3; 
%     vm = 0.3*normrnd(1,0.2); mum = 1/6; betabm = 10*beta; deltam = 1/2;
%     vh = 0; muh = 0; betamh = beta; deltah = 0.6; gammah = 0.06;
%     
%     Nb = Sb + Eb + Ib + Rb;
%     Nm = Sm + Em + Im;
%     
%     Sb1= vb*Nb - (mub + betamb*Im)*Sb + Sb; 
%     Eb1= betamb*Im*Sb - (mub+deltab)*Eb + Eb;
%     Ib1= deltab*Eb - (mub+gammab)*Ib + Ib; 
%     Rb1= gammab*Ib - mub*Rb + Rb;
%     
%     Sm1= vm*Nm -(mum + betabm*Ib)*Sm +Sm;  
%     Em1= betabm*Ib*Sm - (mum + deltam)*Em + Em;
%     Im1= deltam*Em - mum*Im +Im;
%     
%     Sh1= vh - (muh + betamh*Im)*Sh + Sh;
%     Eh1= betamh*Im*Sh - (deltah + gammah + muh)*Eh + Eh;
%     Rh1= deltah*Eh - muh*Rh +Rh;
%     Ch1= gammah*Eh - muh*Ch +Ch;
%     
%     P = [Sb1; Eb1; Ib1; Rb1; Sm1; Em1; Im1; Sh1; Eh1; Rh1; Ch1; beta*normrnd(1,0.2)];
% end
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