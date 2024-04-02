clear; clc; close all
load Days
load beta
load Temp
load Precip

figure 
yyaxis left
plot(Days, beta_t)
ylabel('infection rate')
xlabel('Date')

yyaxis right
plot(Days, T)
ylabel('Temperature')

figure 
yyaxis left
plot(Days, beta_t)
ylabel('infection rate')
xlabel('Date')

yyaxis right
plot(Days, P)
ylabel('Precipitation')