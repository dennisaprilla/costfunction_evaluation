clear; close all; clc;

noise = 3;
trial = 2;
load(sprintf('noise%d_%d.mat', noise, trial));

[X,Y] = meshgrid(r_z, t_z);

figure;
surf(X,Y, rmse_all);
axis equal;
xlabel('Rz (radian)');
ylabel('tz (mm)');
zlabel('RMSE');