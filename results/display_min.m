clear; close all;

filename='allcf_10_1100.mat';
load(filename);

rmse = 1;
gmm  = 2;
gmm2 = 3;

middle = floor(length(r_z)/2);
costfunctions_min_normalized = costfunctions_min - middle;
costfunctions_min_magnitude  = sqrt(sum((costfunctions_min_normalized.^2),2));

figure1 = figure('Name', 'Min Cost Function', 'Position', [50 350 1250 350]);
for i=1:3
    subplot(1,3,i);
    plot(costfunctions_min_normalized(:,1,i), costfunctions_min_normalized(:,2,i), '.r');
    grid on; axis equal;
    xlim([-middle middle]);
    ylim([-middle middle]);
end

figure2 = figure('Name', 'Distance from Zero', 'Position', [20 350 1250 350]);
for i=1:3
    subplot(1,3,i);
    histogram(costfunctions_min_magnitude(:,:,i), 25);
end