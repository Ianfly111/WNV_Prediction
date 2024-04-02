clear; clc; close all
%load Case_GX;

Case_NE = [1,3,7,13,22,39,50,67,84,104,118,130,142];
Pop_NE = 1000000;

Case_CO_ori = [5,31,31,101,102,243,368,438,542,542,584];
Case_CO = Case_CO_ori;
for i = 2:(length(Case_CO)-1)
    Case_CO(i) = 0.5*Case_CO_ori(i) + 0.25*(Case_CO_ori(i-1) + Case_CO_ori(i+1));
end


Pop_CO = 3000000;


y = Case_CO(1:11)/Pop_CO;
Nens = 200;

d = 12;
C = [0 0 0 0 0 0 0 0 0 0 1 0];%; 0 0 0 0 1 1 1 0 0 0 0
q = [1e-4 1e-5 1e-5 1e-7 1e-4 1e-5 1e-5 1e-5 1e-5 1e-5 1e-5 1e-4];
x0 = zeros(d,Nens);
for i=1:Nens
    Eh = normrnd(0.00003,0.00001);
    Sh = 1-Eh;
    Rh = 0;
    Ch = 0;
    
    Sm = normrnd(0.01,0.0001);
    Em = normrnd(0.0001,0.00001);
    Im = normrnd(0.0001,0.00001); %normrnd(0,0.0000001);
    
    Sb = normrnd(0.001,0.0002);
    Eb = normrnd(0.0002,0.00002);
    Ib = normrnd(0.0001,0.00001);
    Rb = normrnd(0.000001,0.000001);
    
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

Nit = 100; 
Est =zeros(d,Nt);
Prt =cell(1,m);
for k = 1:m
    Prt{k} = zeros(d,Nt);
end
Prt_d = zeros(Nit,m);

for nn = 1:Nit
    for ii = 1:Nt
        V = diag(y(:,ii))/(2^ii);% with noise 
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
       Y = y(:,ii)+E;

       %filtering
       K(:,:,ii) = LP_pri*LP_pri'*C'*inv(C*(LP_pri*LP_pri')*C'+Ve); %Gain
       for j = 1:Nens
           X(:,j,ii+1) = max(0,X_pri(:,j,ii) + K(:,:,ii)*(Y(:,j)-C*X_pri(:,j,ii)));
       end
       M(:,ii) = max(0,M_pri(:,ii) + K(:,:,ii)*(y(:,ii)-C*M_pri(:,ii)));
       LP = 1/sqrt(Nens-1)*(X(:,:,ii+1) - M(:,ii));
       P(:,:,ii) = LP*LP';

       Fcst{1}(:,ii) = Forcast(M(:,ii));% +K(:,:,ii)*(y(:,ii)-C*M_pri(:,ii))
       if m>1
           for k = 2:m
               Fcst{k}(:,ii) = Forcast(Fcst{k-1}(:,ii));%+K(:,:,ii)*(y(:,ii)-C*M_pri(:,ii))
           end
       end 
       
    end
    
    Est = Est+M/Nit;
    for k = 1:m
        Prt{k} = Prt{k}+Fcst{k}/Nit;
    end
    
    for k = 1:m
        Prt_d(nn,k) = Fcst{k}(11,ii);
    end
end
%parameter estimation 
figure
subplot(4,1,1)
plot(1:Nt,Est(1,:));
legend('estimated Sb')
xlabel 'Time (weeks)';
ylabel 'Sb'

subplot(4,1,2)
plot(1:Nt,Est(2,:));
legend('estimated Eb')
xlabel 'Time (weeks)';
ylabel 'Eb'

subplot(4,1,3)
plot(1:Nt,Est(3,:));
legend('estimated Ib')
xlabel 'Time (weeks)';
ylabel 'Ib'

subplot(4,1,4)
plot(1:Nt,Est(4,:));
legend('estimated Rb')
xlabel 'Time (weeks)';
ylabel 'Rb'

figure
subplot(3,1,1)
plot(1:Nt,Est(5,:));
legend('estimated Sm')
xlabel 'Time (weeks)';
ylabel 'Sm'

subplot(3,1,2)
plot(1:Nt,Est(6,:));
legend('estimated Em')
xlabel 'Time (weeks)';
ylabel 'Em'

subplot(3,1,3)
plot(1:Nt,Est(7,:));
legend('estimated Im')
xlabel 'Time (weeks)';
ylabel 'Im'


figure
subplot(4,1,1)
plot(1:Nt,Est(8,:));
legend('estimated S_{h}')
xlabel 'Time (weeks)';
ylabel 'Sh'

subplot(4,1,2)
plot(1:Nt,Est(9,:));
legend('estimated E_{h}')
xlabel 'Time (weeks)';
ylabel 'Eh'

subplot(4,1,3)
plot(1:Nt,Est(10,:));
legend('estimated R_{h}')
xlabel 'Time (weeks)';
ylabel 'Rh'

subplot(4,1,4)
plot(1:Nt,Est(11,:));
legend('estimated C_{h}')
xlabel 'Time (weeks)';
ylabel 'Ch'

figure
plot(1:Nt,Est(12,:));
legend('estimated beta')
xlabel 'Time (weeks)';
ylabel 'infection rate'

figure
boxplot(round(Prt_d*Pop_CO),'Symbol',' ');
xlabel 'Time (weeks)';
ylabel 'Distribution'


prediction = round([Prt{1}(11,end), Prt{2}(11,end), Prt{3}(11,end)]*Pop_CO)
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
    
    
    vb = 0.05*normrnd(1,0.2); mub = 0.05; betamb = 5*beta; deltab = 2/3; gammab = 1/3; 
    vm = 0.3*normrnd(1,0.2); mum = 1/6; betabm = 10*beta; deltam = 1/2;
    vh = 0; muh = 0; betamh = beta; deltah = 0.6; gammah = 0.06;
    
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