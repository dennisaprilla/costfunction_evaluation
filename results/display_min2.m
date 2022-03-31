% NOTE: 
% filenaming must be: <bone>_<costfunction>_<param>.mat
% if it is a testset: <bone>_<costfunction>_<param>_testset.mat
% This script is for A-mode simulation. Shifting constant correction is
% implemented in this script
clear; close all;
addpath('..\functions\display\subaxis');

%%

% load the data;
filepath             = 'amode_simulations\accessible_sim2';
filenames            = {'tibia_gmm_scale20', 'tibia_gmm_scale40'};
current_filename     = filenames{2};
use_shiftingconstant = true;
load(strcat(filepath, filesep, current_filename, '.mat'));

% get the information from the trials
noises            = trialsdesc.noises;
pointcounts       = trialsdesc.pointcounts;

% calculate the magnitude (maybe it will be used, or not)
middle = ceil(length(r_z)/2);
costfunctions_min_normalized = costfunctions_min - middle;
costfunctions_min_magnitude  = sqrt(sum((costfunctions_min_normalized.^2),2));

% costfunctions_min contains index of the search-space matrix, let's
% convert it to real rz and tz value (r_z and t_z is available in mat file)
rz_tz_est  = cat(2, r_z(costfunctions_min(:,1,:,:)), t_z(costfunctions_min(:,2,:,:)) );

% if you have shifting constant data, let use_shiftingconstant value to true
if (use_shiftingconstant)
    % load the shifting constant
    newStr = split(current_filename, '_');
    if(length(newStr)==3)
        filename_shiftingconstant = sprintf('%s_shiftingconstant.mat', current_filename);
    else
        temp = strjoin(newStr(1:3), '_');
        filename_shiftingconstant = sprintf('%s_shiftingconstant.mat', temp);
    end
    fullpath_shiftingconstant = strcat(filepath, filesep, filename_shiftingconstant);
    load(fullpath_shiftingconstant);
    % create different name for figure which uses shifting constant
    filename_forsaving = sprintf('%s+sc', current_filename);
else
    filename_forsaving = current_filename;
end

%%

% setup for the figure
figure1       = figure('Name', 'Min Cost Function', 'Position', [50 50 1250 700]);
% figure1       = figure('Name', 'Min Cost Function', 'Position', [50 50 900 700]);
% figure1       = figure('Name', 'Min Cost Function', 'Position', [50 50 250 700]);
subplot_idx   = 1;

% loop for all noises (figure's row)
for noise=1:length(noises)
    
    % loop for all point set configuration (figure's column)
    for pointcount=1:length(pointcounts)
        
        current_pointcount = pointcounts(pointcount);
        
        % is this simulation use shifting constant?
        if (use_shiftingconstant)
            current_shiftingconstant = shiftingconstant(current_pointcount);
        else
            current_shiftingconstant = [0 0]';
        end
        
        % setting up the figure
        subaxis(3,5, subplot_idx, 'Spacing', 0.05, 'MarginLeft',0.075,'MarginRight',0.01,'MarginTop',0.05,'MarginBottom',0.05 );
        % subaxis(3,1, subplot_idx, 'Spacing', 0.05, 'MarginLeft',0.275,'MarginRight',0.040,'MarginTop',0.05,'MarginBottom',0.05 );
        
        % make a scatter plot
        scatter_size = 10;
        rz_est = rz_tz_est(:,1,noise, pointcount) - current_shiftingconstant(1);
        tz_est = rz_tz_est(:,2,noise, pointcount) - current_shiftingconstant(2);
        scatter( rz_est, tz_est, ...    
                 scatter_size, costfunctions_min_magnitude(:, 1, noise, pointcount), 'filled');
        grid on; hold on;
        % draw a rectangle to show the error limit
        rectangle('Position', [-1 -0.001 2 0.002], 'EdgeColor', 'g');
        rectangle('Position', [-2 -0.002 4 0.004], 'EdgeColor', 'r');
        
        % title will be shown only in the first row of figure
        if (noise==1)
            title(sprintf('%d Points', pointcounts(pointcount)));
        end
        
        % ylabel will be shown only in the first column of figure
        if (pointcount==1)
        % if (pointcount==4)
            ylabel(sprintf('Noise %d', noises(noise)), 'fontweight','bold');
        end
        
        % limit the plots so everything is in the same scale
        xlim([r_z(1) r_z(end)]);
        ylim([t_z(1) t_z(end)]);
        
        % increment the subplot index
        subplot_idx = subplot_idx+1;
        
    end
end

% save picture to png
saveas(figure1, sprintf('%s%s%s%s%s', 'pictures', filesep, filepath, filesep, filename_forsaving), 'png');
% save picture to pdf
% https://www.mathworks.com/matlabcentral/answers/12987-how-to-save-a-matlab-graphic-in-a-right-size-pdf
set(figure1,'Units','Inches');
pos = get(figure1,'Position');
set(figure1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(figure1, sprintf('%s%s%s%s%s', 'pictures', filesep, filepath, filesep, filename_forsaving),'-dpdf','-r0');


