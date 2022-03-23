clear; close all;
addpath('..\functions\display\subaxis');

% i forgot to save this from the simulation, so i just put this here
ptCloud_scale = 1000;
range = 15;
step  = 0.5;
r_z   = (-range:step:range);
t_z   = (-range/ptCloud_scale:step/ptCloud_scale:range/ptCloud_scale);

clear ptCloud_scale range step;

%%

% load the data
filepath       = 'abmode_sim1';
fileinfos  = dir(fullfile(filepath, '*.mat'));
filenames  = {fileinfos.name}';

statistical_results = zeros(length(filenames), 6);

for filenumber = 1:length(filenames)
    
    % get the current mat file
    current_filename = filenames{filenumber};
    load(strcat(filepath, filesep, current_filename));
    
    % calculate the magnitude (maybe it will be used, or not)
    % middle = ceil(length(r_z)/2);
    % costfunctions_min_normalized = costfunctions_min - middle;
    % costfunctions_min_magnitude  = sqrt(sum((costfunctions_min_normalized.^2),2));
    
    % costfunctions_min contains index of the search-space matrix, let's
    % convert it to real rz and tz value
    tz_scale   = 1000;
    rz_tz_est  = cat(2, r_z(costfunctions_min(:,1))', t_z(costfunctions_min(:,2))'*tz_scale );
    
    % calculate mean, absolute mean, and variance
    rz_tz_mean    = mean(rz_tz_est, 1);
    rz_tz_meanabs = mean(abs(rz_tz_est), 1);
    rz_tz_cov     = diag(cov(rz_tz_est))';
    
    statistical_results(filenumber, :) = [rz_tz_mean, rz_tz_meanabs, rz_tz_cov];

end

clear costfunctions_min simulation_config fileinfos filenumber current_filename middle ...
      costfunctions_min_normalized costfunctions_min_magnitude tz_scale ...
      rz_tz_est rz_tz_mean rz_tz_meanabs rz_tz_cov;

%% Best Mean

dist_to_origin = vecnorm(statistical_results(:,1:2), 2, 2);
[min_mean_val, min_mean_idx] = min(dist_to_origin);

%% Best Mean Absolute

[min_absmean_rz_val, min_absmean_rz_idx] = min(statistical_results(:,3));
[min_absmean_tz_val, min_absmean_tz_idx] = min(statistical_results(:,4));

%% Best Variance

[min_var_rz_val, min_var_rz_idx] = min(statistical_results(:,5));
[min_var_tz_val, min_var_tz_idx] = min(statistical_results(:,6));

%% Picking the best from our preference
% for example, i will pick the best mean absolute distance for r

% select file
fileidx  = min_var_tz_idx;
filename = filenames{fileidx};
load(strcat(filepath, filesep, filename));

% display the parameter configuration
fprintf('scale_a : %d\nscale_b : %d\nalpha   : %.2f\n', ...
        simulation_config.scale_a, simulation_config.scale_b, simulation_config.alpha);

% calculate the global-minimums position magnitude from the origin
middle = ceil(length(r_z)/2);
costfunctions_min_normalized = costfunctions_min - middle;
costfunctions_min_magnitude  = sqrt(sum((costfunctions_min_normalized.^2),2));

% prepare for scatter plot
scatter_size = 10;
rz_est = r_z(costfunctions_min(:,1))';
tz_est = t_z(costfunctions_min(:,2))';
% draw scatters for global minimums
figure1 = figure('Name', 'Global Minimums', 'Position', [50 50 400 350]);
scatter( rz_est, tz_est, scatter_size, costfunctions_min_magnitude, 'filled');
grid on; hold on;
% draw a rectangle to show the error limit
rectangle('Position', [-1 -0.001 2 0.002], 'EdgeColor', 'g');
rectangle('Position', [-2 -0.002 4 0.004], 'EdgeColor', 'r');
% limit the plots so everything is in the same scale
xlim([r_z(1) r_z(end)]);
ylim([t_z(1) t_z(end)]);
% label the axis
xlabel('R_z');
ylabel('t_z');

% % draw 3d histogram
% figure(2);
% rz_tz_est = [rz_est, tz_est];
% hist3(rz_tz_est, 'Nbins', [10 10], 'CDataMode','auto', 'FaceColor','interp', 'FaceAlpha', 0.8);
% xlabel('R_z');
% ylabel('t_z');
% zlabel('Count');










