clc; close all; clear;
load("processedData.mat");
N=50;
Tplate = zeros(length(parallelTemp)-100,100,N);
TH = zeros(length(parallelTemp)-100,100,N);
TC = zeros(length(parallelTemp)-100,100,N);
TinH = parallelTemp(:,2);
TinC = parallelTemp(:,4);
ToutC = parallelTemp(:,5);
ToutH = parallelTemp(:,3);
for i = 1:N
    TH(1,:,i) = TinH(100)+i*(ToutH(100)-TinH(100))/N;
    TC(1,:,i) = TinC(100)+i*(ToutC(100)-TinC(100))/N;
end
Tplate(1,:,:) = TC(1,:,:);
for i = 1:length(parallelTemp)-101
    [TH(i+1,1,:),Tplate(i+1,1,:),TC(i+1,1,:)] = newTemp(TH(i,100,:),Tplate(i,100,:),TC(i,100,:),TinH(i+100),ToutH(+100),TinC(i+100),ToutC(i+100));
    for j = 1:99
        [TH(i+1,j+1,:),Tplate(i+1,j+1,:),TC(i+1,j+1,:)] = newTemp(TH(i,j,:),Tplate(i,j,:),TC(i,j,:),TinH(i+100),ToutH(i+100),TinC(i+100),ToutC(i+100));
    end
    
end
figure;
dt=1;
time = linspace(100, length(parallelTemp)/dt, length(parallelTemp)-100);
hold on;
% plot(Tplate(:,1,1));
plot(time, TC(:,1,1));
plot(time, TH(:,1,1));
plot(parallelTemp, 'g');
% hold off;
% figure;
% hold on;
% plot(Tplate(:,1,N));
plot(time, TC(:,1,N));
plot(time, TH(:,1,N));
function [newTH,newTplate,newTC] = newTemp(TH,Tplate,TC,TinH,ToutH,TinC,ToutC)

Atot = 0.12; % Total area, m²
Vtot = 0.18e-3; % Total volume, m³
dp = 0.26e-3; % Thickness of conductive panel, m
df = 2.1e-3; % Distance between panels
l = 0.191; % length of panel, m
w = 0.073; % width of panel, m
wTot = 10*w;
dt = 1; % time step, s
N = 50;
alpha = 150; % Convection coefficient, guesstimate
rho = 1e3; % Density of water, kg/m³
rhocopper = 8960; % Density of copper, kg/m³
cp = 4.2e3; % Water cp, J/kg/K
v = 8.9247/2500*1e-3/wTot/df;
sigma = l/v/N; % time constant
lambda = 4e2;
mc = Vtot*rho/2;
mh = Vtot*rho/2;
mp = wTot*l*dp*rhocopper;
cpc = 390; % Copper specific heat, J/kg/K
Apart = Atot/N;
mcpart = mc/N;
mhpart = mh/N;
mppart = mp/N;
% newTH = TH +(-alpha*Atot/mh/cp*(TH-Tplate)+(TinH+ToutH-2*TH)/sigma)*dt/100;
% newTplate = Tplate+ (alpha*Atot/mp/cpc*(TH-2*Tplate+TC))*dt/100;
% newTC = TC + (-alpha*Atot/mc/cp*(TC-Tplate)+mc*cp/sigma*(TinC+ToutC-2*TC))*dt/100;
newTH(1) = TH(1) +(-alpha*Apart/mhpart/cp*(TH(1)-Tplate(1))+(TinH-TH(1))/sigma)*dt/100;
newTH(N) = TH(N) +(-alpha*Apart/mhpart/cp*(TH(N)-Tplate(N))+(ToutH-TH(N))/sigma)*dt/100;
newTC(1) = TC(1) +(-alpha*Apart/mcpart/cp*(TC(1)-Tplate(1))+(TinC-TC(1))/sigma)*dt/100;
newTC(N) = TC(N) +(-alpha*Apart/mcpart/cp*(TC(N)-Tplate(N))+(ToutC-TC(N))/sigma)*dt/100;
newTplate(1) = Tplate(1) + (alpha*Apart/mppart/cpc*(TH(1)-2*Tplate(1)+TC(1))- lambda*dp*wTot/l*N*(Tplate(1)-Tplate(2)))*dt/100;
newTplate(N) = Tplate(N) + (alpha*Apart/mppart/cpc*(TH(N)-2*Tplate(N)+TC(N))- lambda*dp*wTot/l*N*(Tplate(N)-Tplate(N-1)))*dt/100;
for i = 2:N-1
    newTH(i) = TH(i) + (-alpha*Apart/mhpart/cp*(TH(i)-Tplate(i))+(TH(i-1)+TH(i+1)-2*TH(i))/sigma)*dt/100;
    newTC(i) = TC(i) + (-alpha*Apart/mcpart/cp*(TC(i)-Tplate(i))+(TC(i-1)+TC(i+1)-2*TC(i))/sigma)*dt/100;
    newTplate(i) = Tplate(i) + (alpha*Apart/mppart/cpc*(TH(i)-2*Tplate(i)+TC(i))- lambda*dp*wTot/l*N*(2*Tplate(i)-Tplate(i-1)-Tplate(i+1)))*dt/100;
end
end