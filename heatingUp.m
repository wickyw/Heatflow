clc; close all; clear;
data = load("data.mat");
heatingUpData = data.heatinguprenske(:, [1]);
clear data;
%% model 1

totalTime = length(heatingUpData);
heatingCap = 800; %watt
amountOfWater = 10; %liter
density = 0.997;% kg/L
massWater = density*amountOfWater;
time = linspace(0, totalTime, 1);
watt = 2*heatingCap;
specificHeat = 4186; %J/kg grades celsius
heatWater = specificHeat*massWater;
tempWater = zeros(1, totalTime);
tempWater(1:164) = heatingUpData(152);

for i=165:totalTime
    a = watt/heatWater;
    tempWater(i) = tempWater(i-1)+a;
end
%% model 2
tempDiff = heatingUpData(161) - heatingUpData(1);
loss = heatWater*tempDiff/161;
tempWater2 = zeros(1, totalTime);
tempWater2(1:164) = heatingUpData(152);

for i=165:totalTime
    a = (watt+loss)/heatWater;
    tempWater2(i) = tempWater2(i-1)+a;
end

%% plotting
plot(heatingUpData);
xlabel('time [s]')
ylabel(['temperature [', char(176), 'C]']);
title('Water heating up in the water tank, T1');

figure;
plot(heatingUpData);
hold on;
plot(tempWater);
plot(tempWater2);
xlabel('time [s]')
ylabel(['temperature [', char(176), 'C]']);
legend('measurement data', 'first model');
