clc; clear; close all;

addpath(genpath('..\point_cloud_registration_evaluation'));
addpath(genpath("D:\Documents\BELANDA\PhD Thesis\Code\cpp\gmmreg\MATLAB"));

%%
downsample = false;

% read the point cloud (bone) from STL/PLY file
ptCloud          = stlread('data/bone/CT_Tibia_R.stl');
ptCloud_scale    = 1000;
if(downsample)
    ptCloud          = pcdownsample( pointCloud(ptCloud.Points), 'gridAverage', 0.00075 );
    ptCloud_Npoints  = size(ptCloud.Location,1);
    ptCloud_centroid = mean(ptCloud.Location, 1);
    % prepare Ŭ, the noiseless, complete, moving dataset
    U_breve          = (ptCloud.Location - ptCloud_centroid)';
else
    ptCloud_Npoints  = size(ptCloud.Points,1);
    ptCloud_centroid = mean(ptCloud.Points, 1);
    % prepare Ŭ, the noiseless, complete, moving dataset
    U_breve          = (ptCloud.Points - ptCloud_centroid)';
end

% Read the simulated a-mode measurement point cloud, which is a subset of Ŭ.
% These a-mode simulated measurement is manually selected from the bone model.
selectedpoint_str = sprintf('data/bone/amode_measure3.mat');
load(selectedpoint_str);
U = [vertcat(amode_prereg.Position); vertcat(amode_mid.Position)]';
% add isotropic zero-mean gaussian noise to U, simulating noise measuremen
noise        = 1;
% random_noise = normrnd(0, noise/ptCloud_scale, [3, size(U, 2)]);
random_noise = -noise/ptCloud_scale + (noise/ptCloud_scale + noise/ptCloud_scale)*rand(3,size(U, 2));
U            = U + random_noise;

% plot Ŭ, the noiseless, complete, moving dataset
figure1 = figure(1);
figure1.WindowState  = 'maximized';
axes1 = axes('Parent', figure1);
plot3( axes1, ...
       U_breve(1,:), ...
       U_breve(2,:), ...
       U_breve(3,:), ...
       '.', 'Color', [0.7 0.7 0.7],...
       'MarkerSize', 0.1, ...
       'Tag', 'plot_Ubreve');
xlabel('X'); ylabel('Y');
grid(axes1, 'on'); axis(axes1, 'equal'); hold(axes1, 'on');
% plot U, the noisy, incomplete, moving dataset
plot3( axes1, ...
       U(1,:), ...
       U(2,:), ...
       U(3,:), ...
       'or', ...
       'Tag', 'plot_U');
   
%%

model_ptCloud = U';

% obtain all combination of z rotation and translation
range = 5;
step  = 0.25;
r_z   = (-range:step:range);
t_z   = (-range/ptCloud_scale:step/ptCloud_scale:range/ptCloud_scale);
% change z rotation to rotation matrix
rs = [ r_z', zeros(length(r_z), 2) ];
Rs = eul2rotm(deg2rad(rs), 'ZYX');
% change z translation to translation vector
ts = [ zeros(2, length(t_z)); t_z];

% prepare variable to contains all rmse
cf_gmm  = zeros(length(r_z), length(t_z));
cf_gmm2  = zeros(length(r_z), length(t_z));
cf_rmse = zeros(length(r_z), length(t_z));

tic;
for current_z = 1:length(r_z)
    for current_t = 1:length(t_z)
        fprintf('%d %d\n', current_z, current_t);
        
        % transform Ŭ with respected transformation
        U_breve_prime = Rs(:,:,current_z) * U_breve + ts(:,current_t);
        scene_ptCloud = U_breve_prime';
        
%         % RMSE
%         [nearest_idx, nearest_dist] = knnsearch(scene_ptCloud, model_ptCloud);
%         cf_rmse(current_z, current_t) = mean(nearest_dist);
        
        % GMM L2 Distance
        scale = 10e-4;
        [f,~] =  GaussTransform(double(model_ptCloud), double(scene_ptCloud), scale);
        cf_gmm(current_z, current_t) = mean(nearest_dist);

%         % GMM CPD Distance (?)
%         X = scene_ptCloud;
%         Y = model_ptCloud;
%         sigma = 5e-8;
%         w     = 5e-1;
%         [ ~, ~, ~, negativeLogLikelihood ] = computeEStep(X, Y, sigma, w);
        cf_gmm2(current_z, current_t) = mean(nearest_dist);

        
%         % display what is happening
%         delete(findobj('Tag', 'plot_Ubreveprime'));
%         plot3( axes1, ...
%                U_breve_prime(1,:), ...
%                U_breve_prime(2,:), ...
%                U_breve_prime(3,:), ...
%                '.g', 'MarkerSize', 0.1, ...
%                'Tag', 'plot_Ubreveprime');
%         drawnow;
        
    end
end
toc;

%%
[X,Y] = meshgrid(r_z, t_z);

figure2 = figure(2);
figure2.WindowState  = 'maximized';

subplot(1,2,1);
surf(X,Y, cf_rmse);
xlabel('Rz (deg)');
ylabel('tz (mm)');
zlabel('GMM L2 distance');
view(90, 90);

subplot(1,2,2);
surf(X,Y, cf_gmm2);
xlabel('Rz (deg)');
ylabel('tz (mm)');
zlabel('GMM CPD distance');
view(90, 90);

