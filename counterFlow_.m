clc; close all; clear;
load("processedData.mat");
N=99;
dt=1;
counterTemp_ = counterTemp;
Tplate = zeros(length(counterTemp)-100/dt,100,N);
TH = zeros(length(counterTemp)-100/dt,100,N);
TC = zeros(length(counterTemp)-100/dt,100,N);
TinH = counterTemp(:,2);
TinC = counterTemp(:,5);
ToutC = counterTemp(:,4);
ToutH = counterTemp(:,3);
Ttank(1) = parallelTemp(100, 1);
amountWater = 0.01; %L/s
WaterDt = amountWater*dt;
WaterRho = 0.997; %kg/L
massWaterDt = WaterDt*WaterRho;
for i = 1:N
    TH(1,:,i) = TinH(100)+i*(ToutH(100)-TinH(100))/N;
    TC(1,:,N+1-i) = TinC(100)+i*(ToutC(100)-TinC(100))/N;
end
Tplate(1,:,:) = TC(1,:,:);
for i = 1:length(counterTemp)/dt-101
    [TH(i+1,1,:),Tplate(i+1,1,:),TC(i+1,1,:)] = newTemp(TH(i,100,:),Tplate(i,100,:),TC(i,100,:),TinH(i+100),ToutH(i+100),TinC(i+100),ToutC(i+100));
    for j = 1:99
        [TH(i+1,j+1,:),Tplate(i+1,j+1,:),TC(i+1,j+1,:)] = newTemp(TH(i,j,:),Tplate(i,j,:),TC(i,j,:),TinH(i+100),ToutH(i+100),TinC(i+100),ToutC(i+100));
    end
    Ttank(i+1) = (Ttank(i)*WaterRho*(10-WaterDt)+(TH(i+1, 100, N)*WaterRho*WaterDt))/(WaterRho*10);
end
%dt = 1;
time = linspace(100, length(counterTemp)/dt, length(counterTemp)-100);
figure;
hold on;
% plot(Tplate(:,1,1));
plot(time, TC(:,1,1), 'b');
plot(time, TH(:,1,1), 'k');
plot(counterTemp_, 'g');
% hold off;
% figure;
% hold on;
% plot(Tplate(:,1,N));
plot(time, TC(:,1,N), 'b');
plot(time, TH(:,1,N), 'k');
plot(time,Ttank);
function [newTH,newTplate,newTC] = newTemp(TH,Tplate,TC,TinH,ToutH,TinC,ToutC)

Atot = 0.12; % Total area, m²
Vtot = 0.18e-3; % Total volume, m³
dp = 0.26e-3; % Thickness of conductive panel, m
df = 2.1e-3; % Distance between panels
l = 0.191; % length of panel, m
w = 0.073; % width of panel, m
wTot = 10*w;
dt = 1; % time step, s
N = 99;
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
newTC(1) = TC(1) +(-alpha*Apart/mcpart/cp*(TC(1)-Tplate(1))+(ToutC-TC(1))/sigma)*dt/100;
newTC(N) = TC(N) +(-alpha*Apart/mcpart/cp*(TC(N)-Tplate(N))+(TinC-TC(N))/sigma)*dt/100;
newTplate(1) = Tplate(1) + (alpha*Apart/mppart/cpc*(TH(1)-2*Tplate(1)+TC(1))- lambda*dp*wTot/l*N*(Tplate(1)-Tplate(2)))*dt/100;
newTplate(N) = Tplate(N) + (alpha*Apart/mppart/cpc*(TH(N)-2*Tplate(N)+TC(N))- lambda*dp*wTot/l*N*(Tplate(N)-Tplate(N-1)))*dt/100;
for k = 2:N-2
    newTH(k) = TH(k) + (-alpha*Apart/mhpart/cp*(TH(k)-Tplate(k))+(TH(k-1)-TH(k))/sigma)*dt/100;
    newTC(N-k) = TC(N-k) + (-alpha*Apart/mcpart/cp*(TC(N-k)-Tplate(N-k))+(TC(N-k+1)-TC(N-k))/sigma)*dt/100;
    newTplate(k) = Tplate(k) + (alpha*Apart/mppart/cpc*(TH(k)-2*Tplate(k)+TC(k))- lambda*dp*wTot/l*N*(2*Tplate(k)-Tplate(k-1)-Tplate(k+1)))*dt/100;
end
end