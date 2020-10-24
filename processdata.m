clc; close all; clear;
load('data.mat');
plot(counterheat(:,[1,3]), counterheat(:,[2,4]));
figure;
plot(countertemp(:,1:2:9), countertemp(:,2:2:10));