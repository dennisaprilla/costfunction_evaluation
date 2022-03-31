clear; close all;
addpath('..\functions\display\subaxis');

% % i forgot to save this from the simulation, so i just put this here
% ptCloud_scale = 1000;
% range = 15;
% step  = 0.5;
% r_z   = (-range:step:range);
% t_z   = (-range/ptCloud_scale:step/ptCloud_scale:range/ptCloud_scale);
% clear ptCloud_scale range step;

savefile = false;

%% Loading the File

% load the data
simname     = 'abmode_sim2d';
simpath     = 'abmode_simulations';
simfullpath = strcat(simpath, filesep, simname);
fileinfos   = dir(fullfile(simfullpath, '*.mat'));
filenames   = {fileinfos.name}';

statistical_results = zeros(length(filenames), 6);

for filenumber = 1:length(filenames)
    
    % get the current mat file
    current_filename = filenames{filenumber};
    load(strcat(simfullpath, filesep, current_filename));
    
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

%% Summarizing the results

% Best Mean
dist_to_origin = vecnorm(statistical_results(:,1:2), 2, 2);
[min_val, min_idx] = min(dist_to_origin);
statistical_results_preference(1).name     = 'mean_rztz';
statistical_results_preference(1).fullname = 'Mean Rz-tz';
statistical_results_preference(1).idx      = min_idx;
statistical_results_preference(1).val      = min_val;

% Best Mean Absolute
[min_val, min_idx] = min(statistical_results(:,3));
statistical_results_preference(2).name     = 'absmean_rz';
statistical_results_preference(2).fullname = 'Absolute Mean Rz';
statistical_results_preference(2).idx      = min_idx;
statistical_results_preference(2).val      = min_val;
[min_val, min_idx] = min(statistical_results(:,4));
statistical_results_preference(3).name     = 'absmean_tz';
statistical_results_preference(3).fullname = 'Absolute Mean tz';
statistical_results_preference(3).idx      = min_idx;
statistical_results_preference(3).val      = min_val;

% Best Variance
[min_val, min_idx] = min(statistical_results(:,5));
statistical_results_preference(4).name     = 'var_rz';
statistical_results_preference(4).fullname = 'Absolute Var Rz';
statistical_results_preference(4).idx      = min_idx;
statistical_results_preference(4).val      = min_val;
[min_val, min_idx] = min(statistical_results(:,6));
statistical_results_preference(5).name     = 'var_tz';
statistical_results_preference(5).fullname = 'Absolute Var tz';
statistical_results_preference(5).idx      = min_idx;
statistical_results_preference(5).val      = min_val;

clear min_val min_idx dist_to_origin;

%% Picking the best from our preference
% for example, i will pick the best mean absolute distance for rz

MEAN_RZTZ  = 1;
ABSMEAN_RZ = 2;
ABSMEAN_TZ = 3;
VAR_RZ     = 4;
VAR_TZ     = 5;

% select file
preference = VAR_RZ;
fileidx    = statistical_results_preference(preference).idx;
filename   = filenames{fileidx};
load(strcat(simfullpath, filesep, filename));

% display the parameter configuration
fprintf('scale_a : %d\nscale_b : %d\nalpha   : %.2f\n', ...
        simulation_config.scale_a, simulation_config.scale_b, simulation_config.alpha);
    
% display the statistical results based on preference
fprintf('mean_rz\t\t: %.4f\tmean_tz\t\t: %.4f\nabsmean_rz\t: %.4f\tabsmean_tz\t: %.4f\nvar_rz\t\t: %.4f\tvar_tz\t\t: %4f\n', ...
        statistical_results(fileidx,1), statistical_results(fileidx,2), ...
        statistical_results(fileidx,3), statistical_results(fileidx,4), ...
        statistical_results(fileidx,5), statistical_results(fileidx,6) );
        
% for copying to excel
final_result = [ statistical_results(fileidx,3), ...
                 statistical_results(fileidx,4), ...
                 statistical_results(fileidx,1), ...
                 statistical_results(fileidx,2), ...
                 statistical_results(fileidx,5), ...
                 statistical_results(fileidx,6), ...
               ];
    
%% Displaying the Result

%{
figure1 = figure('Name', 'Global Minimums', 'Position', [50 50 375 350]);
% prepare for scatter plot
scatter_size = 10;
rz_est = r_z(costfunctions_min(:,1))';
tz_est = t_z(costfunctions_min(:,2))';
% draw scatters for global minimums
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
%}

figure1 = figure('Name', 'Global Minimums', 'Position', [50 50 500 650]);
% prepare histogram
rz_est = r_z(costfunctions_min(:,1))';
tz_est = t_z(costfunctions_min(:,2))';
% draw 3d histogram
hist3( [rz_est, tz_est], ...
       'Ctrs', {-5:0.25:5, -5/1000:0.25/1000:5/1000}, ...
       'FaceAlpha', 0.5, ...
       'EdgeColor', 'interp', ...
       'CDataMode','auto', 'FaceColor','interp');
hold on;

% i want to put intensity map under the 3d histogram, so i follow this
% example https://www.mathworks.com/help/stats/hist3.html at "Plot
% Histogram with Intensity Map"
maxbin_z            = 8;
zlevel_intensitymap = -8;

% count the number of each element in each bin
N = hist3([rz_est, tz_est], 'Ctrs', {r_z, t_z}, 'CDataMode','auto', 'FaceColor','interp', 'FaceAlpha', 0.8 );
% Generate a grid to draw the 2-D projected view of intensities by using pcolor.
N_pcolor = N';
N_pcolor(size(N_pcolor,1)+1,size(N_pcolor,2)+1) = 0;
xl = linspace(min(r_z),max(r_z),size(N_pcolor,2));
yl = linspace(min(t_z),max(t_z),size(N_pcolor,1));
% Draw the intensity map by using pcolor. 
h = pcolor(xl, yl, N_pcolor);
colormap('default'); 
colorbar;
caxis([0 maxbin_z])
% Set the z-level of the intensity map to view the histogram and the intensity map together.
% zlevel_intensitymap = -max(N_pcolor(:));
h.ZData = zlevel_intensitymap * ones(size(N_pcolor));
h.EdgeColor = 'none';
ax = gca;
ax.ZTick(ax.ZTick < 0) = [];

% draw the rectangle to show the accuracy limit
plot3([-1 -1 1 1 -1], [-0.001 0.001 0.001 -0.001 -0.001], [zlevel_intensitymap zlevel_intensitymap zlevel_intensitymap zlevel_intensitymap zlevel_intensitymap], 'MarkerFaceColor', 'g');
plot3([-2 -2 2 2 -2], [-0.002 0.002 0.002 -0.002 -0.002], [zlevel_intensitymap zlevel_intensitymap zlevel_intensitymap zlevel_intensitymap zlevel_intensitymap], 'MarkerFaceColor', 'r');

view(-30, 20);
zlim([zlevel_intensitymap maxbin_z]);
xlabel('R_z');
ylabel('t_z');
zlabel('Count');
title(sprintf('%s %s %s_%d', simname, statistical_results_preference(preference).name, 'noise', simulation_config.noise.level), 'Interpreter', 'none');

clear N N_pcolor xl yl h ax;

%% Save Pictures

if (savefile)

% get directory
[file, path] = uiputfile('*.*', 'File Selection', sprintf('%s_%s%d', statistical_results_preference(preference).name, 'noise', simulation_config.noise.level));
% if user press cancel, stop
if file==0
  return
end

% save picture as png
saveas(figure1, strcat(path, file), 'png');
% save picture to pdf
% https://www.mathworks.com/matlabcentral/answers/12987-how-to-save-a-matlab-graphic-in-a-right-size-pdf
set(figure1,'Units','Inches');
pos = get(figure1,'Position');
set(figure1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(figure1, strcat(path, file),'-dpdf','-r0');
% save 
writestruct(simulation_config, strcat(path, file, '.xml'));

end








