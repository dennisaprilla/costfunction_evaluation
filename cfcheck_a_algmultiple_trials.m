clc; clear; close all;

addpath(genpath('..\pointcloudregistration_evaluations'));
addpath(genpath('..\gmmreg\MATLAB'));

displaybone = true;
displaycf   = true;

%%

% read the point cloud (bone) from STL/PLY file
% ptCloud          = stlread('data/bone/CT_Femur_R.stl');
ptCloud          = stlread('data/bone/CT_Tibia_R.stl');
ptCloud_scale    = 1000;
ptCloud_Npoints  = size(ptCloud.Points,1);
ptCloud_centroid = mean(ptCloud.Points, 1);
% prepare Ŭ, the noiseless, complete, moving dataset
U_breve          = (ptCloud.Points - ptCloud_centroid)';

% (for debugging only) show figure for sanity check
if(displaybone)
    figure1 = figure('Name', 'Bone', 'Position', [0 -100 400 900]);
    axes1 = axes('Parent', figure1);
    plot3( axes1, ...
           U_breve(1,:), ...
           U_breve(2,:), ...
           U_breve(3,:), ...
           '.', 'Color', [0.7 0.7 0.7], ...
           'MarkerSize', 0.1, ...
           'Tag', 'plot_bone_full');
    xlabel('X'); ylabel('Y'); zlabel('Z');
end

% Read the simulated a-mode measurement point cloud, which is a subset of Ŭ.
% These a-mode simulated measurement is manually selected from the bone model.
selectedpoint_str = sprintf('data/bone/amode_accessible_sim2/amode_tibia_20.mat');
load(selectedpoint_str);
U = vertcat(amode_all.Position)';

% (for debugging only) show figure for sanity check
if(displaybone)
    grid on; axis equal; hold on;
    plot3( axes1, ...
           U(1,:), ...
           U(2,:), ...
           U(3,:), ...
           'or', 'MarkerFaceColor', 'r', ...
           'Tag', 'plot_bone_full');
end
    
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

% trials configuration
num_trials        = 100;
noise             = 2;
num_costfunction  = 3;
costfunctions_min = ones(num_trials, 2, num_costfunction); 

for trial=1:num_trials
        
    fprintf('trials: %d\n', trial);
    
    % add isotropic zero-mean gaussian noise to U, simulating noise measuremen
    % random_noise = normrnd(0, noise/ptCloud_scale, [3, size(U, 2)]);
    random_noise  = -noise/ptCloud_scale + (noise/ptCloud_scale + noise/ptCloud_scale)*rand(3,size(U, 2));
    model_ptCloud = (U + random_noise)';
    
    % (for debugging only) show figure for sanity check
    if (displaybone)
        grid on; axis equal; hold on;
        plot3( axes1, ...
               model_ptCloud(:,1), ...
               model_ptCloud(:,2), ...
               model_ptCloud(:,3), ...
               'oy', ...
               'Tag', 'plot_bone_full');
    end
       
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
            U_breve_prime = Rs(:,:,current_z) * U_breve + ts(:,current_t);
            scene_ptCloud = U_breve_prime';
            
            % RMSE
            [nearest_idx, nearest_dist] = knnsearch(scene_ptCloud, model_ptCloud);
            cf_t_rmse(current_t) = mean(nearest_dist);
            
            % GMM L2 Distance
            scale = 40e-4;
            [f,~] =  GaussTransform(double(model_ptCloud), double(scene_ptCloud), scale);
            cf_t_gmm(current_t) = -f;
    
            % GMM CPD Distance (?)
            X = scene_ptCloud;
            Y = model_ptCloud;
            % sigma = 5e-8; w = 5e-1;
            % sigma = 9e-8; w = 1e-90;
            sigma = 1e-7; w = 1e-10;
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
    
    % (for debugging only) display cost function surface
    if (displaycf)
        [X,Y] = meshgrid(r_z, t_z);
        figure2 = figure(2);
        figure2.WindowState  = 'maximized';
        subplot(1,3,1);
        surf(X,Y, cf_rmse);
        view([-10, 10]);
        title('RMSE');
        xlabel('Rz (deg)');
        ylabel('tz (mm)');
        zlabel('Cost');
        subplot(1,3,2);
        surf(X,Y, cf_gmm);
        view([-10, 10]);
        title(sprintf('GMM L2 Distance'));
        xlabel('Rz (deg)');
        ylabel('tz (mm)');
        zlabel('Cost');
        subplot(1,3,3);
        surf(X,Y, cf_gmm2);
        view([-10, 10]);
        title(sprintf('GMM Loglikelihood'));
        xlabel('Rz (deg)');
        ylabel('tz (mm)');
        zlabel('Cost');
    end
    
    % (for debugging only) break the loop
    if (or(displaybone, displaycf))
        break;
    end

end

% save('results\rmsegmm_tibia30_2.mat', 'costfunctions_min', 'r_z', 't_z');
