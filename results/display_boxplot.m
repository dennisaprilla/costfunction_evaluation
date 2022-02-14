clear;
addpath(genpath('..\functions\display'));

% filename='tibia_gmm_scale20';
% filename='tibia_gmm_scale30';
filename='tibia_gmm_scale40';
% filename='tibia_rmse';
load(strcat(filename, '.mat'));

% setup the noises (this is supposed to be given from the mat file, but i
% forgot to save it.
noises            = [1 2 3];
pointcounts       = [15 20 25 30];

% calculate the magnitude (maybe it will be used, or not)
middle = ceil(length(r_z)/2);
costfunctions_min_normalized = costfunctions_min - middle;
costfunctions_min_magnitude  = sqrt(sum((costfunctions_min_normalized.^2),2));

% costfunctions_min contains index of the search-space matrix, let's
% convert it to real rz and tz value
rz_tz_est  = cat(2, r_z(costfunctions_min(:,1,:,:)), t_z(costfunctions_min(:,2,:,:))*1000 );
rz_tz_abs  = abs(rz_tz_est);
rz_tz_mean = mean(rz_tz_abs, 1);
rz_tz_std  = std(rz_tz_abs, 1);
disp(rz_tz_std(:,:,:,1));

% % rearrange data, the requirement for boxplotGroup, please refer
% % https://www.mathworks.com/matlabcentral/answers/331381-how-do-you-create-a-grouped-boxplot-with-categorical-variables-on-the-x-axis#answer_418952
noise_asPrimary = 1;
data={};

if(noise_asPrimary==1)
    for dof=1:2
        for current_noise=1:length(noises)
            data{dof, current_noise} = reshape( rz_tz_abs(:, dof, current_noise, :), [], length(pointcounts));
        end
    end
else
    for dof=1:2
        for curent_pointcounts=1:length(pointcounts)
            data{dof, curent_pointcounts} = reshape( rz_tz_abs(:, dof, :, curent_pointcounts), [], length(noises));
        end
    end
end

% we use subaxis function to control more for the spacing for the subplot
figure1 = figure('Name', 'Translation Error', 'Position', [0 0 1200 450]);
for dof=1:2
    
    subplot(1,2, dof); hold on;
    
    data_temp = data(dof, :);
    if (noise_asPrimary==1)
        primarylabel   = {'1', '2', '3'};
        secondarylabel = {'15', '20', '25', '30'};
    else
        primarylabel   = {'15', '20', '25', '30'};
        secondarylabel = {'1', '2', '3'};
    end
    
    boxplot_handle{dof} = boxplotGroup( data_temp, 'PrimaryLabels', primarylabel, 'SecondaryLabels', secondarylabel); 
    grid on;
    
    if(dof==1)
        ylabel('Absolute Difference (degree)');
        title('Rotation Error');
    else
        ylabel('Absolute Difference (mm)');
        title('Translation Error');
    end
    
    ylim([0, 10]);
end












