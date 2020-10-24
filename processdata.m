clc; close all; clear;
% load('data.mat');
% plot(counterheat(:,[1,3]), counterheat(:,[2,4]));
% figure;
% plot(countertemp(:,1:2:9), countertemp(:,2:2:10));

clc; close all; clear;
data = load("data.mat");
counterFlow = data.counterflow(:,[2,4]);
counterHeat = data.counterheat(:,[2 4]);
counterTemp = data.countertemp(:,[2 4 6 8 10]);
eersteMetingFlow = data.eerstemetingflow(:,[2,4]);
eersteMetingHeat = data.eerstemetingheat(:,[2,4]);
eersteMetingTemp = data.eerstemetingtemp(:,[2 4 6 8 10]);
heatingupRenskeNew = data.eerstemetingtemp(:,[2 4 6 8 10]); 
parallelFlow = data.parallelflow(:,[2 4]);
parallelHeat = data.parallelheat(:,[2 4]);
parallelTemp = data.paralleltemp(:,[2 4 6 8 10]);
clear data;
save("processedData.mat"); 