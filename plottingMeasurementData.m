clear all; close all; clc;

data = load("processedData.mat");
countertemp_ = data.counterTemp;
paralleltemp_ = data.parallelHeat;
plot(countertemp_);
xlabel('time [s]');
ylabel(['temperature [', char(176), 'C]']);
title('counter flow');
figure;
plot(paralleltemp_);
xlabel('time [s]');
ylabel(['temperature [', char(176), 'C]']);
title('parallel flow');


clear data;