clear; close all;

filename='allcf_amode3_3.mat';
load(filename);

rmse = 1;
gmm  = 2;
gmm2 = 3;

middle = ceil(length(r_z)/2);
costfunctions_min_normalized = costfunctions_min - middle;
costfunctions_min_magnitude  = sqrt(sum((costfunctions_min_normalized.^2),2));

rz_tz_est  = cat(2, r_z(costfunctions_min(:,1,:)), t_z(costfunctions_min(:,2,:)) );
rz_tz_mean = mean(abs(rz_tz_est), 1);
rz_tz_median = median(abs(rz_tz_est), 1);

%%

figure1 = figure('Name', 'Min Cost Function', 'Position', [50 350 1250 350]);
for i=1:3
    subplot(1,3,i);
%     plot(costfunctions_min_normalized(:,1,i), costfunctions_min_normalized(:,2,i), '.r');
%     grid on; axis equal;
%     xlim([-middle middle]);
%     ylim([-middle middle]);
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

figure2 = figure('Name', 'Distance from Zero', 'Position', [70 350 1250 350]);
for i=1:3
    subplot(1,3,i);
    histogram(costfunctions_min_magnitude(:,:,i), 10);
end

%%
figure3 = figure('Name', '3D Histogram', 'Position', [90 350 1250 350]);
for i=1:3
    subplot(1,3,i);
    hist3(rz_tz_est(:,:,i), 'Nbins', [10 10], 'CDataMode','auto', 'FaceColor','interp');
end