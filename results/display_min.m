clear; close all;

% filename='allcf_amode4_1_inlier.mat';
filename='rmsegmm_tibia30_2.mat';
load(filename);

rmse = 1;
gmm  = 2;
gmm2 = 3;

middle = ceil(length(r_z)/2);
costfunctions_min_normalized = costfunctions_min - middle;
costfunctions_min_magnitude  = sqrt(sum((costfunctions_min_normalized.^2),2));

rz_tz_est  = cat(2, r_z(costfunctions_min(:,1,:)), t_z(costfunctions_min(:,2,:)) );
rz_tz_mean = mean(abs(rz_tz_est), 1);
rz_tz_std = std(abs(rz_tz_est), 1);
disp("mean");
disp(rz_tz_mean);
disp("std");
disp(rz_tz_std);

%%

figure1 = figure('Name', 'Min Cost Function', 'Position', [20 100 1500 440]);
for i=1:3

    subplot(1,3,i);
    plot(rz_tz_est(:,1,i), rz_tz_est(:,2,i), '.b');
    grid on; hold on;

    rectangle('Position', [-1 -0.001 2 0.002], 'EdgeColor', 'g')
    rectangle('Position', [-2 -0.002 4 0.004], 'EdgeColor', 'r')

    xlim([r_z(1) r_z(end)]);
    ylim([t_z(1) t_z(end)]);
    xlabel('Rotation Error (deg)');
    ylabel('Translation Error (mm)');
    
end

%%

% figure2 = figure('Name', 'Distance from Zero', 'Position', [70 350 1250 350]);
% for i=1:3
%     subplot(1,3,i);
%     histogram(costfunctions_min_magnitude(:,:,i), 10);
% end

%%
figure3 = figure('Name', '3D Histogram', 'Position', [20 100 1500 500]);
for i=1:3
    subplot(1,3,i);
    hist3(rz_tz_est(:,:,i), 'Nbins', [10 10], 'CDataMode','auto', 'FaceColor','interp', 'FaceAlpha', 0.8);
end
