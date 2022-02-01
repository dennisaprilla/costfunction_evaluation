clc; clear; close all;

addpath(genpath('..\pointcloudregistration_evaluations'));
addpath(genpath('..\gmmreg\MATLAB'));

%% Prepare Bone Point Cloud and A-mode Measurement

% read the point cloud (bone) from STL/PLY file
ptCloud          = stlread('data/bone/CT_Femur_R.stl');
% ptCloud          = stlread('data/bone/CT_Tibia_R.stl');
ptCloud_scale    = 1000;
ptCloud_Npoints  = size(ptCloud.Points,1);
ptCloud_centroid = mean(ptCloud.Points, 1);
% prepare Ŭ, the noiseless, complete, moving dataset
U_breve          = (ptCloud.Points - ptCloud_centroid)';

% Read the simulated a-mode measurement point cloud, which is a subset of Ŭ.
% These a-mode simulated measurement is manually selected from the bone model.
selectedpoint_str = sprintf('data/bone/amode_measure4.mat');
load(selectedpoint_str);
U = [vertcat(amode_prereg.Position); vertcat(amode_mid.Position)]';

clear ptCloud;

%% Prepare Inlier Assumption

threshold = 20 / ptCloud_scale;
[nearest_idx, nearest_dist] = knnsearch(U_breve', U', 'k', 1000);

U_breve_inliers = [];
for i=1:size(U,2)
    U_breve_inliers = [U_breve_inliers, U_breve(:, nearest_idx(i,nearest_dist(i,:)<threshold))];
end
U_breve_inliers = unique(U_breve_inliers', 'rows')';

clear U_breve;

%% Simulate Search Space
    
% obtain all combination of z rotation and translation
range = 10;
step  = 0.25;
r_z   = (-range:step:range);
t_z   = (-range/ptCloud_scale:step/ptCloud_scale:range/ptCloud_scale);
% change z rotation to rotation matrix
rs = [ r_z', zeros(length(r_z), 2) ];
Rs = eul2rotm(deg2rad(rs), 'ZYX');
% change z translation to trgitanslation vector
ts = [ zeros(2, length(t_z)); t_z];

%% Trials

num_trials        = 100;
noise             = 1;
num_costfunction  = 3;
costfunctions_min = zeros(num_trials, 2, num_costfunction); 

for trial=1:num_trials
        
    fprintf('trials: %d\n', trial);
    
    % add isotropic zero-mean gaussian noise to U, simulating noise measuremen
    % random_noise = normrnd(0, noise/ptCloud_scale, [3, size(U, 2)]);
    random_noise  = -noise/ptCloud_scale + (noise/ptCloud_scale + noise/ptCloud_scale)*rand(3,size(U, 2));
    model_ptCloud = (U + random_noise)';
       
    %%
    
    % prepare variable to contains all rmse
    cf_gmm  = zeros(length(r_z), length(t_z));
    cf_gmm2 = zeros(length(r_z), length(t_z));
    cf_rmse = zeros(length(r_z), length(t_z));
    
    loop_z = length(r_z);
    loop_t = length(t_z);
    
    tic;
    parfor current_z = 1:loop_z
        
        cf_t_rmse = zeros(1, loop_t);
        cf_t_gmm  = zeros(1, loop_t);
        cf_t_gmm2 = zeros(1, loop_t);
        
        for current_t = 1:loop_t
            
            % transform Ŭ with respected transformation
            % U_breve_prime = Rs(:,:,current_z) * U_breve + ts(:,current_t);
            % scene_ptCloud = U_breve_prime';
            U_breve_prime_inliers = Rs(:,:,current_z) * U_breve_inliers + ts(:,current_t);
            scene_ptCloud = U_breve_prime_inliers';
            
            % RMSE
            [nearest_idx, nearest_dist] = knnsearch(scene_ptCloud, model_ptCloud);
            cf_t_rmse(current_t) = mean(nearest_dist);
            
            % GMM L2 Distance
%             scale = 40e-4;
            scale = 37.5e-4;
            [f,~] =  GaussTransform(double(model_ptCloud), double(scene_ptCloud), scale);
            cf_t_gmm(current_t) = -f;
    
            % GMM CPD Distance (?)
            X = scene_ptCloud;
            Y = model_ptCloud;
%             sigma = 5e-8;
%             w     = 5e-1;
%             sigma = 9e-8;
%             w     = 1e-90;
%             sigma = 1e-7;
%             w     = 1e-10;
            sigma = 5e-7;
            w     = 0.9;
            [ ~, ~, ~, negativeLogLikelihood ] = computeEStep(X, Y, sigma, w);
            cf_t_gmm2(current_t) = negativeLogLikelihood;
            
            
        end
        
        cf_rmse(current_z, : ) = cf_t_rmse;
        cf_gmm(current_z, : )  = cf_t_gmm;
        cf_gmm2(current_z, : ) = cf_t_gmm2;
    end
    
    toc;
    
    minValue = min(cf_rmse(:));
    [costfunctions_min(trial, 1, 1), costfunctions_min(trial, 2, 1)] = find(cf_rmse == minValue);
    minValue = min(cf_gmm(:));
    [costfunctions_min(trial, 1, 2), costfunctions_min(trial, 2, 2)] = find(cf_gmm == minValue);
    minValue = min(cf_gmm2(:));
    [costfunctions_min(trial, 1, 3), costfunctions_min(trial, 2, 3)] = find(cf_gmm2 == minValue);
    
    %%
    
%     [X,Y] = meshgrid(r_z, t_z);
% 
%     figure2 = figure(2);
%     figure2.WindowState  = 'maximized';
%     subplot(1,3,1);
%     surf(X,Y, cf_rmse);
%     title('RMSE');
%     xlabel('Rz (deg)');
%     ylabel('tz (mm)');
%     zlabel('Cost');
%     subplot(1,3,2);
%     surf(X,Y, cf_gmm);
%     title(sprintf('GMM L2 Distance'));
%     xlabel('Rz (deg)');
%     ylabel('tz (mm)');
%     zlabel('Cost');
%     subplot(1,3,3);
%     surf(X,Y, cf_gmm2);
%     title(sprintf('GMM Loglikelihood'));
%     xlabel('Rz (deg)');
%     ylabel('tz (mm)');
%     zlabel('Cost');

end

save('results\allcf_amode4_1_inlier.mat', 'costfunctions_min', 'r_z', 't_z');
