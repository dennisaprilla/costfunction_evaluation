load('tibia_gmm_scale10_600.mat');
costfunctions_min1 = costfunctions_min;
load('tibia_gmm_scale10_150.mat');
costfunctions_min2 = costfunctions_min;

costfunctions_min = cat(1, costfunctions_min1, costfunctions_min2);
save('tibia_gmm_scale10.mat', 'costfunctions_min', 'r_z', 't_z');

% load('tibia_rmse_500.mat');
% costfunctions_min1 = costfunctions_min;
% load('tibia_rmse_250.mat');
% costfunctions_min2 = costfunctions_min;
% 
% costfunctions_min = cat(1, costfunctions_min1, costfunctions_min2);
% save('tibia_rmse.mat', 'costfunctions_min', 'r_z', 't_z');