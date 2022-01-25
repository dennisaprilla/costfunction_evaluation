clear; close all;

filename='test.mat';
load(filename);


rmse = 1;
gmm  = 2;
gmm2 = 3;

figure;
for i=1:3
    subplot(1,3,i);
    scatter(costfunctions_min(:,1,i), costfunctions_min(:,2,i));
    xlim([0 41]);
    ylim([0 41]);
    grid on; axis equal;
end