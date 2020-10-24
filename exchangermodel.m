clc; close all; clear;
%% Model of the Heat exchanger
% This model is based on the theory described in the MOD04 and MOD05
% reader, chapter 14 and section 16.4.2.
load("processedData.mat");
A = 0.012; % Area per plate, m²
Atot = 0.12; % Total area, m²
V = 0.018e-3; % Volume per plate, m³
Vtot = 0.18e03; % Total volume, m³
dp = 0.26e-3; % Thickness of conductive panel, m
df = 2.1e-3; % Distance between panels
l = 0.191; % length of panel, m
w = 0.073; % width of panel, m
np = 10; % number of plates
N = 10; % number of model segments per plate
dl = l/N; % length of model segment
dt = 5; % time step, s
tmax = length(parallelTemp)*dt; % end time, s
v = 1; % Fluid flow rate, m/s
sigma = dl/v; % time constant
cflow = dt/sigma; % constant for convective heat transport
alpha = 25; % Convection coefficient
lambda = 400;
Rth = (dp/lambda+2/alpha)/A*N;
rho = 1e3; % Density of water, kg/m³
mWaterSegment = V*rho;
cp = 4.2e3; % Water cp, J/kg/K

%% Model data arrays initialisation:
Tmodel = zeros(tmax/dt,np+1,N);
Tin = zeros(length(parallelTemp),2*np-1);
for i = 1:2:np+1
    Tin(:,i) = parallelTemp(:,3);
    Tmodel(1,:,1:N/2) = parallelTemp(1,3);
    Tmodel(1,:,(N/2+1):N) = parallelTemp(1,2);
end
for i = 2:2:np+1
    Tin(:,i) = parallelTemp(:,5);
    Tmodel(1,:,1:N/2) = parallelTemp(1,5);
    Tmodel(1,:,(N/2+1):N) = parallelTemp(1,4);
end
dTflow = zeros(1,np+1,N);
dTplate = zeros(1,np+1,N); % Temperature change through interaction through plate
for i = 1:tmax/dt-1
    dTplate(1,1,1:N) = -(Tmodel(i,1,1:N)-Tmodel(i,2,1:N))/Rth*dt;
    for j = 2:np
        % First, convective heat transport through the fluid:
        dTflow(1,j,1) = cflow*(Tin(i,j)-Tmodel(i,j,1));
        dTflow(1,j,2:N) = cflow*(Tmodel(i,j,1:N-1)-Tmodel(i,j,2:N));
        % Second, through-plate heat transport
        dTplate(1,j,1:N) = -(2*Tmodel(i,j,1:N)-Tmodel(i,j+1,1:N)-Tmodel(i,j-1,1:N))/Rth*dt;
    end
    dTplate(1,np+1,N) = -(Tmodel(i,j,N)-Tmodel(i,j-1,N))/Rth*dt/mWaterSegment/cp;
    dT = dTflow + dTplate;
    Tmodel(i+1,:,:) = Tmodel(i,:,:) + dT;
end