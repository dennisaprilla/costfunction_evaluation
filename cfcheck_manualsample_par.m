clc; clear; close all;

addpath(genpath('..\point_cloud_registration_evaluation'));
addpath(genpath("D:\Documents\BELANDA\PhD Thesis\Code\cpp\gmmreg\MATLAB"));

%%

% read the point cloud (bone) from STL/PLY file
ptCloud          = stlread('data/bone/CT_Tibia_R.stl');
ptCloud_scale    = 1000;
ptCloud_Npoints  = size(ptCloud.Points,1);
ptCloud_centroid = mean(ptCloud.Points, 1);
% prepare Ŭ, the noiseless, complete, moving dataset
U_breve          = (ptCloud.Points - ptCloud_centroid)';

% Read the simulated a-mode measurement point cloud, which is a subset of Ŭ.
% These a-mode simulated measurement is manually selected from the bone model.
selectedpoint_str = sprintf('data/bone/amode_measure.mat');
load(selectedpoint_str);
U = [vertcat(amode_prereg.Position); vertcat(amode_mid.Position)]' ./ ptCloud_scale;
% add isotropic zero-mean gaussian noise to U, simulating noise measuremen
noise        = 2;
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
% arrange view 
view(axes1, 5, 80);
   
%%

model_ptCloud = U';

% obtain all combination of z rotation and translation
range = 10;
step  = 0.5;
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

loop_z = length(r_z);
loop_t = length(t_z);

tic;
parfor current_z = 1:loop_z
    
    cf_t_rmse = zeros(1, loop_t);
    cf_t_gmm2 = zeros(1, loop_t);
    
    for current_t = 1:loop_t
        
%         fprintf('%d %d\n', current_z, current_t);
        
        % transform Ŭ with respected transformation
        U_breve_prime = Rs(:,:,current_z) * U_breve + ts(:,current_t);
        scene_ptCloud = U_breve_prime';
        
        % compute nearest index (and nearest distance) using knnsearch
        [nearest_idx, nearest_dist] = knnsearch(scene_ptCloud, model_ptCloud);
        % store the mean distance
        cf_t_rmse(current_t) = mean(nearest_dist);
        
%         scale = 0.00075;
%         [f,~] =  GaussTransform(double(model_ptCloud), double(scene_ptCloud), scale);
%         cf_gmm(current_z, current_t) = -f;

        X = scene_ptCloud;
        Y = model_ptCloud;
        sigma = 1e-7;
        w     = 1e-8;
        [ ~, ~, ~, negativeLogLikelihood ] = computeEStep(X, Y, sigma, w);
        cf_t_gmm2(current_t) = negativeLogLikelihood;
        
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
    
    cf_rmse(current_z, : ) = cf_t_rmse;
    cf_gmm2(current_z, : ) = cf_t_gmm2;
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

