clear; close all

addpath('..\functions\display\subaxis');

% filename='tibia_gmm_scale20';
% filename='tibia_gmm_scale30';
% filename='tibia_gmm_scale40';
filename='tibia_rmse';
load(strcat(filename, '.mat'));

noises            = [1 2 3];
pointcounts       = [15 20 25 30];

middle = ceil(length(r_z)/2);
costfunctions_min_normalized = costfunctions_min - middle;
costfunctions_min_magnitude  = sqrt(sum((costfunctions_min_normalized.^2),2));

rz_tz_est  = cat(2, r_z(costfunctions_min(:,1,:,:)), t_z(costfunctions_min(:,2,:,:)) );

n_noises    = size(costfunctions_min, 3);
n_numpoints = size(costfunctions_min, 4);
figure1     = figure('Name', 'Min Cost Function', 'Position', [50 50 900 700]);
subplot_idx = 1;

for current_noise=1:n_noises
    for current_numpoints=1:n_numpoints
        subaxis(3,4, subplot_idx, 'Spacing', 0.05, 'MarginLeft',0.075,'MarginRight',0.01,'MarginTop',0.05,'MarginBottom',0.05 );
        
        scatter_size = 10;
        scatter( rz_tz_est(:,1,current_noise, current_numpoints), rz_tz_est(:,2,current_noise, current_numpoints), ...
                 scatter_size, costfunctions_min_magnitude(:, 1, current_noise, current_numpoints), 'filled');
        grid on; hold on;
        rectangle('Position', [-1 -0.001 2 0.002], 'EdgeColor', 'g');
        rectangle('Position', [-2 -0.002 4 0.004], 'EdgeColor', 'r');
        
        if (current_noise==1)
            title(sprintf('%d Points', pointcounts(current_numpoints)));
        end
        
        if (current_numpoints==1)
            ylabel(sprintf('Noise %d', noises(current_noise)), 'fontweight','bold');
        end
        
        xlim([r_z(1) r_z(end)]);
        ylim([t_z(1) t_z(end)]);
        
        subplot_idx = subplot_idx+1;
        
    end
end

% save picture to png
saveas(figure1, sprintf('pictures/%s', filename), 'png');
% save picture to pdf
% https://www.mathworks.com/matlabcentral/answers/12987-how-to-save-a-matlab-graphic-in-a-right-size-pdf
set(figure1,'Units','Inches');
pos = get(figure1,'Position');
set(figure1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(figure1, sprintf('pictures/%s', filename),'-dpdf','-r0');