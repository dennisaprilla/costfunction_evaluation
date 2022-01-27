clc; clear; close all;

addpath(genpath('..\point_cloud_registration_evaluation'));
addpath(genpath("D:\Documents\BELANDA\PhD Thesis\Code\cpp\gmmreg\MATLAB"));

%%

% ptCloud            = pcread('data/bunny/reconstruction/bun_zipper.ply');
% ptCloud            = pcdownsample( ptCloud, 'gridAverage', 0.01);
% ptCloud_scale      = 1000;
% ptCloud_Npoints    = size(ptCloud.Location, 1);
% ptCloud_centroid   = mean(ptCloud.Location, 1);
% U_breve            = (ptCloud.Location - ptCloud_centroid)';
% 
% U             = ( pcdownsample( pointCloud(U_breve'), 'random', 0.03).Location )';
% noise         = 3;
% random_noise  = normrnd(0, noise/ptCloud_scale, [3, size(U, 2)]);
% U             = U + random_noise;

% read the point cloud (bone) from STL/PLY file
ptCloud          = stlread('data/bone/CT_Tibia_R.stl');
ptCloud_scale    = 1000;
ptCloud_Npoints  = size(ptCloud.Points,1);
ptCloud_centroid = mean(ptCloud.Points, 1);
% prepare Ŭ, the noiseless, complete, moving dataset
U_breve          = (ptCloud.Points - ptCloud_centroid)';

%{
% sampling point from preregistration area, this while iteration is
% making sure the selected sample points in each preregistration area
selectedpoint_str = 'data/bone/amode_areasphere_443cm.mat';
load(selectedpoint_str);
% is not too near to each other
preregistrationArea_number = length(preregistrationArea);
preregistrationArea_samplingNumber = 4;
accepted_condition = false;
threshold_distance = 10;
while(~accepted_condition)
    preregistrationArea_samplePoints = {};
    for i=1:preregistrationArea_number
        point_number = size(preregistrationArea{i}, 1);

        preregistrationArea_randomIndex = randi( [1, point_number], 1, preregistrationArea_samplingNumber );
        preregistrationArea_samplePoints{i} = preregistrationArea{i}(preregistrationArea_randomIndex, :);

        samplePoints_index = 1:preregistrationArea_samplingNumber;
        combination_index = nchoosek(samplePoints_index,2);

        for j=1:length(combination_index)
            points_toCheck = [ preregistrationArea_samplePoints{i}(combination_index(j,1), :); ...
                               preregistrationArea_samplePoints{i}(combination_index(j,2), :) ];

            if ( pdist( points_toCheck, 'euclidean') < (threshold_distance/ptCloud_scale) )
                accepted_condition = false;
                break;
            else
                accepted_condition = true;
            end

        end

        if (accepted_condition == false)
            break;
        end
    end
end

% gather all the amode sample points to single array
amode_samplePoints = [];
for i=1:preregistrationArea_number
    amode_samplePoints = [ amode_samplePoints; preregistrationArea_samplePoints{i} ];
end

% combine the preregistration area from cell to matrix
preregistrationArea = cell2mat(preregistrationArea');
%}


% Read the simulated a-mode measurement point cloud, which is a subset of Ŭ.
% These a-mode simulated measurement is manually selected from the bone model.
% selectedpoint_str = sprintf('data/bone/amode_measure.mat');
% load(selectedpoint_str);
% amode_samplePoints = [ amode_samplePoints; vertcat(amode_mid.Position)/ptCloud_scale];
selectedpoint_str = sprintf('data/bone/amode_measure.mat');
load(selectedpoint_str);
amode_samplePoints = [vertcat(amode_prereg.Position)/ptCloud_scale; vertcat(amode_mid.Position)/ptCloud_scale];

% add noise to scene point cloud
noise              = 0;
% random_noise       = normrnd(0, noise/ptCloud_scale, [size(amode_samplePoints, 1), 3]);
random_noise       = -noise/ptCloud_scale + (noise/ptCloud_scale + noise/ptCloud_scale)*rand(size(amode_samplePoints, 1), 3);
% amode_samplePoints = amode_samplePoints + random_noise;


%%

% % display the femure bone
% figure1 = figure(1);
% figure1.WindowState  = 'maximized';
% axes1 = axes('Parent', figure1);
% plot3( axes1, ...
%        U_breve(:,1), ...
%        U_breve(:,2), ...
%        U_breve(:,3), ...
%        '.', 'Color', [0.7 0.7 0.7], ...
%        'MarkerSize', 0.1, ...
%        'Tag', 'plot_bone_full');
% xlabel('X'); ylabel('Y'); zlabel('Z');
% grid on; axis equal; hold on;
% % plot a-mode sample points
% plot3( axes1, ...
%        amode_samplePoints(:,1), ...
%        amode_samplePoints(:,2), ...
%        amode_samplePoints(:,3), ...
%        'oy', 'MarkerFaceColor', 'y',...
%        'Tag', 'plot_bone_samplepoints');

   
%%

% renaming variable
model_ptCloud = amode_samplePoints;

% obtain all combination of z rotation and translation
range  = 5;
step   = 0.5;
r_z   = (-range:step:range);
t_z   = (-range/ptCloud_scale:step/ptCloud_scale:range/ptCloud_scale);
% change z rotation to rotation matrix
rs = [ r_z', zeros(length(r_z), 2) ];
Rs = eul2rotm(deg2rad(rs), 'ZYX');
% change z translation to translation vector
ts = [ zeros(2, length(t_z)); t_z];

% prepare variable to contains all rmse
cf_gmm  = zeros(length(r_z), length(t_z));
cf_rmse = zeros(length(r_z), length(t_z));

for current_z = 1:length(r_z)
    for current_t = 1:length(t_z)
        fprintf('%d %d\n', current_z, current_t);
       
        % transform Ŭ with respected transformation
        U_breve_prime = Rs(:,:,current_z) * U_breve + ts(:,current_t);
        scene_ptCloud = U_breve_prime';

        % compute nearest index (and nearest distance) using knnsearch
        [nearest_idx, nearest_dist] = knnsearch(scene_ptCloud, model_ptCloud);
        % store the mean distance
        cf_rmse(current_z, current_t) = mean(nearest_dist);  
        
    end
end

%%

[X,Y] = meshgrid(r_z, t_z);

figure2 = figure(2);
figure2.WindowState  = 'maximized';

% subplot(1,2,1);
% surf(X,Y, cf_gmm);
% xlabel('Rz (deg)');
% ylabel('tz (mm)');
% zlabel('GMM L2 distance');

subplot(1,2,2);
surf(X,Y, cf_rmse);
xlabel('Rz (deg)');
ylabel('tz (mm)');
zlabel('RMSE');

minValue = min(cf_gmm(:));
[row1, column1] = find(cf_gmm == minValue);
