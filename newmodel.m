clc; close all; clear;
load("processedData.mat");
Atot = 0.12; % Total area, m²
Vtot = 0.18e03; % Total volume, m³
dp = 0.26e-3; % Thickness of conductive panel, m
df = 2.1e-3; % Distance between panels
l = 0.191; % length of panel, m
w = 0.073; % width of panel, m
wTot = 10*w;
np = 10; % number of plates
dt = 5; % time step, s
tmax = length(parallelTemp)*dt; % end time, s
v = wTot*df; % Fluid flow rate, m/s
sigma = l/v; % time constant
cflow = dt/sigma; % constant for flow heat transport
alpha = 25; % Convection coefficient, guesstimate
lambda = 400; % Conduction coefficient of copper
rho = 1e3; % Density of water, kg/m³
cp = 4.2e3; % Water cp, J/kg/K
N = 10;
Tplate = zeros(length(parallelTemp),1);
TH = zeros(length(parallelTemp),1);
TC = zeros(length(parallelTemp),1);
TinH = parallelTemp(:,5);
TinC = parallelTemp(:,3);
ToutC = parallelTemp(:,2);
ToutH = parallelTemp(:,4);
TH(1) = (TinH(1)+ToutH(1))/2;
TC(1) = (TinC(1)+ToutC(1))/2;
Tplate(1) = (TH(1)+TC(1))/2;
for i = 1:length(parallelTemp)-1
    TH(i+1) = TH(i) -alpha*Atot/m/cp*(TH(i)-Tplate(i))+(TinH(i)-TH(i))/sigma-(TH(i)-ToutH(i))/sigma;
    
end