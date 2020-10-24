clc; close all; clear;
%% Model of the Heat exchanger
% This model is based on the theory described in the MOD04 and MOD05
% reader, chapter 14 and section 16.4.2.
A = 0.012; % Area per plate, m²
Atot = 0.12; % Total area, m²
V = 0.018e-3; % Volume per plate, m³
Vtot = 0.18e03; % Total volume, m³
m = 0.65; % Heat exchanger mass, kg
dp = 0.26*10^(-3); % Thickness of conductive panel, m
df = 2.1*10^(-3); % Distance between panels
l = 0.191; % length of panel, m
w = 0.073; % width of panel, m
np = 10; % number of plates
N = 100; % number of model segments per plate
dl = l/N; % length of model segment
dt = 5; % time step, s
tmax = 1000; % end time, s
v = 0; % Fluid flow rate, m/s
labda = 401; %thermal heat conductivity of copper
TempCooling = 20;
amountWater = 0.0036; %L/s
WaterDt = amountWater*dt;
WaterRho = 0.997; %kg/L
massWaterDt = WaterDt*WaterRho;
heatCapWater = 4200;
%% alleen warme stroom, plaat, koude stroom
koudeStroom = zeros(tmax/dt, N);
koudeStroom(:,1) = TempCooling;
warmeStroom = zeros(tmax/dt, N);

warmeStroom(1,:) = 60;
koudeStroom(1,:) = TempCooling;
tempPlaat = 17;
AdeelPlaat = A/N;
Tnew= 60;
for i = 2:tmax/dt
    warmeStroom(i,1) = Tnew;
    for j = 2:N
%         if j == 1
%             warmeStroom(i, j) = 60;
%             koudeStroom(i, j) = TempCooling;
%         else
            a = labda*AdeelPlaat/dp;
%             warmeStroom(i-1, j-1)
%             koudeStroom(i-1, j-1)
            Qwarm = -(labda*AdeelPlaat/dp)*dt*dl*(warmeStroom(i-1, j-1)-koudeStroom(i-1, j-1));
%             Qwarm/(massWaterDt*heatCapWater)
            Qkoud = -Qwarm;
            warmeStroom(i, j) = Qwarm/(massWaterDt*heatCapWater) + warmeStroom(i-1, j-1);
            koudeStroom(i, j) = Qkoud/(massWaterDt*heatCapWater) + koudeStroom(i-1, j-1);
%         end
        
    end
    Tnew = (Tnew*WaterRho*(10-WaterDt)+warmeStroom(i, N)*WaterRho*WaterDt)/(WaterRho*10)
end
%% Model data arrays initialisation:
Tmodel = zeros(tmax/dt,np,N);
for i = 1:tmax/dt
    for j = 1:np
        for k = 1:N
            % First, convective heat transport through the fluid:
            
            
            
            
        end
    end
end