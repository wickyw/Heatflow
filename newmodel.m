clc; close all; clear;
load("processedData.mat");
Atot = 0.12; % Total area, m²
Vtot = 0.18e-3; % Total volume, m³
dp = 0.26e-3; % Thickness of conductive panel, m
df = 2.1e-3; % Distance between panels
l = 0.191; % length of panel, m
w = 0.073; % width of panel, m
wTot = 10*w;
np = 10; % number of plates
dt = 5; % time step, s
tmax = length(parallelTemp)*dt; % end time, s
alpha = 30; % Convection coefficient, guesstimate
lambda = 400; % Conduction coefficient of copper
rho = 1e3; % Density of water, kg/m³
rhocopper = 8960; % Density of copper, kg/m³
cp = 4.2e3; % Water cp, J/kg/K
Tplate = zeros(length(parallelTemp),100);
TH = zeros(length(parallelTemp),100);
TC = zeros(length(parallelTemp),100);
TinH = parallelTemp(:,2);
TinC = parallelTemp(:,4);
ToutC = parallelTemp(:,5);
ToutH = parallelTemp(:,3);
TH(1,:) = (TinH(1)+ToutH(1))/2;
TC(1,:) = (TinC(1)+ToutC(1))/2;
Tplate(1,:) = (TH(1)+TC(1))/2;
v = mean(mean(parallelFlow))/2500*1e-3/wTot/df;
sigma = l/v; % time constant
mc = Vtot*rho/2;
mh = Vtot*rho/2;
mp = wTot*l*dp*rhocopper;
cpc = 390; % Copper specific heat, J/kg/K
for i = 1:length(parallelTemp)-1
    if mod(i,100) == 0
        i;
    end
    TH(i+1,1) = TH(i,100)+(-alpha*Atot/mh/cp*(TH(i,100)-Tplate(i,100))+(TinH(i+1)-TH(i,100))/sigma-(TH(i,100)-ToutH(i+1))/sigma)*dt/100;
    Tplate(i+1,1) = Tplate(i,100) +(alpha*Atot/mp/cpc*(TH(i,100)-2*Tplate(i,100)+TC(i,100)))*dt/100;
    TC(i+1,1) = TC(i,100) +(-alpha*Atot/mc/cp*(TC(i,100)-Tplate(i,100))+mc*cp/sigma*(TinC(i+1)+ToutC(i+1)-2*TC(i,100)))*dt/100;
    for j = 1:99
        TH(i+1,j+1) = TH(i,j)+(-alpha*Atot/mh/cp*(TH(i,j)-Tplate(i,j))+(TinH(i)+(TinH(i+1)-TinH(i))*j/100-TH(i,j))/sigma-(TH(i,j)-ToutH(i)-(ToutH(i+1)-ToutH(i))*j/100)/sigma)*dt/100;
        Tplate(i+1,j+1) = Tplate(i,j) +(alpha*Atot/mp/cpc*(TH(i,j)-2*Tplate(i,j)+TC(i,j)))*dt/100;
        TC(i+1,j+1) = TC(i,j) +(-alpha*Atot/mc/cp*(TC(i,j)-Tplate(i,j))+mc*cp/sigma*(TinC(i)+(TinC(i+1)-TinC(i))*j/100+ToutC(i)+(ToutC(i+1)-ToutC(i))*j/100-2*TC(i,j)))*dt/100;
    end
    
end
figure;
hold on;
plot(Tplate(:,1));
plot(TC(:,1));
plot(TH(:,1));