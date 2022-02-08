% load('tibia_gmm_scale40_500.mat');
% costfunctions_min1 = costfunctions_min;
% load('tibia_gmm_scale40_150.mat');
% costfunctions_min2 = costfunctions_min;
% load('tibia_gmm_scale40_100.mat');
% costfunctions_min3 = costfunctions_min;
% 
% costfunctions_min = cat(1, costfunctions_min1, costfunctions_min2, costfunctions_min3);
% save('tibia_gmm_scale40.mat', 'costfunctions_min', 'r_z', 't_z');

load('tibia_rmse_500.mat');
costfunctions_min1 = costfunctions_min;
load('tibia_rmse_250.mat');
costfunctions_min2 = costfunctions_min;

costfunctions_min = cat(1, costfunctions_min1, costfunctions_min2);
save('tibia_rmse.mat', 'costfunctions_min', 'r_z', 't_z');