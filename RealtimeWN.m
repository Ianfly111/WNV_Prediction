clear; clc; close all

load Case_GX; % pop 2770000 use 10:50 
Pop = 2770000;
CumCase = cumsum(Case_GX(10:45))/Pop;
y = CumCase(1:30);
d = 11;
Nens = 200;
x0 = zeros(d,Nens);
for i=1:Nens
    Sh = normrnd(0.9999,0.00001);
    Eh = normrnd(0.00003,0.00001);
    Ih = normrnd(0.00004,0.00001);
    if Ih < 0
        Ih = 0;
    end
    Rh = 1-Sh-Eh-Ih;
    if Rh < 0
        Rh = 0;
        Ih = 1-Sh;
    end
    
    Sm = normrnd(0.01,0.0001);
    Em = normrnd(0.0001,0.00001);
    Im = normrnd(0.0001,0.00001); %normrnd(0,0.0000001);
    
    Sb = normrnd(0.001,0.0002);
    Eb = normrnd(0.0002,0.00002);
    Ib = normrnd(0.0001,0.00001);
    Rb = normrnd(0.000001,0.000001);
    
    x0(:,i) = [Sb, Eb, Ib, Rb, Sm, Em, Im, Sh, Eh, Ih, Rh];%,BetaM
end
% **********************************************%
C = [0 0 0 0 0 0 0 0 0 0 1];%; 0 0 0 0 1 1 1 0 0 0 0
[No,Nt] = size(y);
X_pri = zeros(d,Nens,Nt);
M_pri = zeros(d,Nt);
X = zeros(d,Nens,Nt+1);
M = zeros(d,Nt);
P_pri = zeros(d,d,Nt);
P = zeros(d,d,Nt);
K = zeros(d,No,Nt);
X(:,:,1) = x0;

m = 5;
Fcst = cell(1,m);  
for k = 1:m
    Fcst{k} = zeros(d,Nt);
end

Nit = 1; 
Est =zeros(d,Nt);
Prt =cell(1,m);
for k = 1:m
    Prt{k} = zeros(d,Nt);
end

for nn = 1:Nit
    for ii = 1:Nt
        V = diag(y(:,ii))/(2^ii);% with noise 
        q = [1e-4 1e-5 1e-5 1e-7 1e-4 1e-5 1e-5 1e-5 1e-5 1e-5 1e-5]/(ii);% with noise
        Q = diag(q);
   %prediction
       for j = 1:Nens
           X_pri(:,j,ii) = Forcast(X(:,j,ii));
       end 

       X_pri(:,:,ii) = X_pri(:,:,ii) + Q*randn(size(X_pri(:,:,ii)));
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
           X(:,j,ii+1) = X_pri(:,j,ii) + K(:,:,ii)*(Y(:,j)-C*X_pri(:,j,ii));
       end
       M(:,ii) = M_pri(:,ii) + K(:,:,ii)*(y(:,ii)-C*M_pri(:,ii));
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
plot(1:Nt,Est(6,:));
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
plot(1:Nt,Est(6,:));
legend('estimated Im')
xlabel 'Time (weeks)';
ylabel 'Im'


figure
subplot(4,1,1)
plot(1:Nt,Est(8,:));
legend('estimated S_{h}')
xlabel 'Time (weeks)';
ylabel 'Suceptible human'

subplot(4,1,2)
plot(1:Nt,Est(9,:));
legend('estimated E_{h}')
xlabel 'Time (weeks)';
ylabel 'Eh'

subplot(4,1,3)
plot(1:Nt,Est(10,:));
legend('estimated I_{h}')
xlabel 'Time (weeks)';
ylabel 'Ih'

subplot(4,1,4)
plot(1:Nt,Est(11,:));
legend('estimated R_{h}')
xlabel 'Time (weeks)';
ylabel 'Removed human'


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
    Sh = SEIRSEISEIR(8); Eh = SEIRSEISEIR(9); Ih = SEIRSEISEIR(10); Rh = SEIRSEISEIR(11);
    
    vb = 0.1; mub = 0.05; betamb = 0.01; deltab = 2/3; gammab = 1/3; 
    vm = 0.3; mum = 1/6; betabm = 0.01; deltam = 1/2;
    vh = 0; muh = 0; betamh = 0.2; deltah = 1; gammah = 1;
    
    Sb1= vb - (mub + betamb*Im)*Sb + Sb; 
    Eb1= betamb*Im*Sb - (mub+deltab)*Eb + Eb;
    Ib1= deltab*Eb - (mub+gammab)*Ib + Ib; 
    Rb1= gammab*Ib - mub*Rb + Rb;
    
    Sm1= vm -(mum + betabm*Ib)*Sm +Sm;  
    Em1= betabm*Ib*Sm - (mum + deltam)*Em + Em;
    Im1= deltam*Em - mum*Im +Im;
    
    Sh1= vh - (muh + betamh*Im)*Sh + Sh;
    Eh1= betamh*Im*Sh - (muh+deltah)*Eh + Eh;
    Ih1= deltah*Eh - (muh + gammah)*Ih + Ih;
    Rh1= gammah*Ih -muh*Rh + Rh;
    
    P = [Sb1; Eb1; Ib1; Rb1; Sm1; Em1; Im1; Sh1; Eh1; Ih1; Rh1];
end