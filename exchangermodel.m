clc; close all; clear;
%% Model of the Heat exchanger
% This model is based on the theory described in the MOD04 and MOD05
% reader, chapter 14 and section 16.4.2.
A = 0.012; % Area per plate, m²
Atot = 0.12; % Total area, m²
V = 0.018e-3; % Volume per plate, m³
Vtot = 0.18e03; % Total volume, m³
dp = 0.26e-3; % Thickness of conductive panel, m
df = 2.1e-3; % Distance between panels
l = 0.191; % length of panel, m
w = 0.073; % width of panel, m
np = 10; % number of plates
N = 100; % number of model segments per plate
dl = l/N; % length of model segment
dt = 5; % time step, s
tmax = 100; % end time, s
v = 5; % Fluid flow rate, m/s
sigma = dl/v; % time constant
cflow = dt/sigma; % constant for convective heat transport
alpha = 10; % Convection coefficient
lambda = 400;
Rth = (dp/lambda+2/alpha)*1/A/N;
rho = 1e3; % Density of water, kg/m³
mWaterSegment = dl*w*df*rho;
cp = 4.2e3; % Water cp, J/kg/K
%% Model data arrays initialisation:
Tmodel = zeros(tmax/dt,np,N);
dTflow = zeros(N,1);
dTplate = zeros(N,1); % Temperature change through interaction through plate
for i = 1:tmax/dt-1
    for j = 1:np        
        % First, convective heat transport through the fluid:
        dTflow(1) = cflow*(Tin(i)-Tmodel(i,j,1));
        dTflow(2:N) = cflow*(Tmodel(i,j,1:N-1)-Tmodel(i,j,2:N));
        % Second, convective heat transport 
        dTplate(1) = -(Tmodel(i,j,1)-Tmodel(i,j+1,1))/Rth*dt;
        dTplate(2:N-1) = -(2*Tmodel(i,j,2:N-1)-Tmodel(i,j+1,2:N-1)-Tmodel(i,j-1,2:N-1))/Rth*dt;
        dTplate(N) = -(Tmodel(i,j,N)-Tmodel(i,j-1,N))/Rth*dt/mWaterSegment/cp;
        Tmodel(i+1,j,:) = dTflow + dTplate;
    end
end