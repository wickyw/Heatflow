clc; close all; clear;
load("processedData.mat");
Tplate = zeros(length(parallelTemp),100);
TH = zeros(length(parallelTemp),100);
TC = zeros(length(parallelTemp),100);
TinH = parallelTemp(:,1);
TinC = parallelTemp(:,4);
ToutC = parallelTemp(:,5);
ToutH = parallelTemp(:,3);
TH(1,:) = TinH(1);
TC(1,:) = (TinC(1)+ToutC(1))/2;
Tplate(1,:) = TC(1);
for i = 1:length(parallelTemp)-1
    [TH(i+1,1),Tplate(i+1,1),TC(i+1,1)] = newTemp(TH(i,100),Tplate(i,100),TC(i,100),TinH(i),ToutH(i),TinC(i),ToutC(i));
    for j = 1:99
        [TH(i+1,j+1),Tplate(i+1,j+1),TC(i+1,j+1)] = newTemp(TH(i,j),Tplate(i,j),TC(i,j),TinH(i),ToutH(i),TinC(i),ToutC(i));
    end
    
end
figure;
hold on;
plot(Tplate(:,1));
plot(TC(:,1));
plot(TH(:,1));
function [newTH,newTplate,newTC] = newTemp(TH,Tplate,TC,TinH,ToutH,TinC,ToutC)

Atot = 0.12; % Total area, m²
Vtot = 0.18e-3; % Total volume, m³
dp = 0.26e-3; % Thickness of conductive panel, m
df = 2.1e-3; % Distance between panels
l = 0.191; % length of panel, m
w = 0.073; % width of panel, m
wTot = 10*w;
dt = 5; % time step, s
N = 10;
alpha = 30; % Convection coefficient, guesstimate
rho = 1e3; % Density of water, kg/m³
rhocopper = 8960; % Density of copper, kg/m³
cp = 4.2e3; % Water cp, J/kg/K
v = 8.9247/2500*1e-3/wTot/df;
sigma = l/v; % time constant
mc = Vtot*rho/2;
mh = Vtot*rho/2;
mp = wTot*l*dp*rhocopper;
cpc = 390; % Copper specific heat, J/kg/K
newTH = TH +(-alpha*Atot/mh/cp*(TH-Tplate)+(TinH+ToutH-2*TH)/sigma)*dt/100;
newTplate = Tplate+ (alpha*Atot/mp/cpc*(TH-2*Tplate+TC))*dt/100;
newTC = TC + (-alpha*Atot/mc/cp*(TC-Tplate)+mc*cp/sigma*(TinC+ToutC-2*TC))*dt/100;

end