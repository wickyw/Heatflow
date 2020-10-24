clc; close all; clear;
%% Model of the Heat exchanger
% This model is based on the theory described in the MOD04 and MOD05
% reader, chapter 14 and section 16.4.2.
A = 0.012; % Area per plate, m²
Atot = 0.12; % Total area, m²
V = 0.018e-3; % Volume per plate, m³
Vtot = 0.18e03; % Total volume, m³
m = 0.65; % Heat exchanger mass, kg
dp = 0.26e-3; % Thickness of conductive panel, m
df = 2.1e-3; % Distance between panels
l = 0.191; % length of panel, m
w = 0.073; % width of panel, m
np = 10; % number of plates
N = 100; % number of model segments per plate
dl = l/N; % length of model segment
dt = 5; % time step, s
tmax = 100; % end time, s
v = 0; % Fluid flow rate, m/s

%% Model data arrays initialisation:
Tmodel = zeros(tmax/dt,np,N);
for i = 1:tmax/dt
    for j = 1:np
        for k = 1:N
            % First, convective heat transport through the fluid:
            
            
            
            
        end
    end
end