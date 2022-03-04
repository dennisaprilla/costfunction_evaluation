clc; clear; close all;

addpath(genpath('..\pointcloudregistration_evaluations'));
addpath(genpath('..\gmmreg\MATLAB'));

displaybone = true;
displaycf   = true;

%%

% read the point cloud (bone) from STL/PLY file
ptCloud          = stlread('data/bone/CT_Tibia_R.stl');
% ptCloud          = stlread('data/bone/CT_Femur_R.stl');
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
step  = 0.5;
r_z   = (-range:step:range);
t_z   = (-range/ptCloud_scale:step/ptCloud_scale:range/ptCloud_scale);
% change z rotation to rotation matrix
rs = [ r_z', zeros(length(r_z), 2) ];
Rs = eul2rotm(deg2rad(rs), 'ZYX');
% change z translation to translation vector
ts = [ zeros(2, length(t_z)); t_z];

% trials configuration
num_trials   = 1;
noise        = 2;
costfunction_name = "gmm1";
costfunctions_min = zeros(num_trials, 2); 

for trial=1:num_trials
    
    fprintf('trials: %d\n', trial);

    % add isotropic zero-mean gaussian noise to U, simulating noise measuremen
    % random_noise = normrnd(0, noise/ptCloud_scale, [3, size(U, 2)]);
    random_noise = -noise/ptCloud_scale + (noise/ptCloud_scale + noise/ptCloud_scale)*rand(3,size(U, 2));
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
    cf  = zeros(length(r_z), length(t_z));

    loop_z = length(r_z);
    loop_t = length(t_z);

    tic;
    parfor current_z = 1:loop_z

        cf_t  = zeros(1, loop_t);

        for current_t = 1:loop_t

            % transform Ŭ with respected transformation
            U_breve_prime = Rs(:,:,current_z) * U_breve + ts(:,current_t);
            scene_ptCloud = U_breve_prime';
            
            if(strcmp(costfunction_name, "gmm1"))
                % GMM L2 Distance
                scale = 20e-4;
                [f,~] =  GaussTransform(double(model_ptCloud), double(scene_ptCloud), scale);
                cf_t(current_t) = -f;
                
            elseif (strcmp(costfunction_name, "gmm2"))
                X = scene_ptCloud;
                Y = model_ptCloud;
                % sigma = 9e-8; w = 1e-90;
                % sigma = 5e-7; w = 0.9;
                sigma = 5e-8; w = 5e-1;
                [ ~, ~, ~, negativeLogLikelihood ] = computeEStep(X, Y, sigma, w);
                cf_t(current_t) = negativeLogLikelihood;

            elseif (strcmp(costfunction_name, "rmse"))
                % RMSE
                [nearest_idx, nearest_dist] = knnsearch(scene_ptCloud, model_ptCloud);
                cf_t(current_t) = mean(nearest_dist);
            end

        end

        cf(current_z, : )  = cf_t;

    end
    toc;

    % (for debugging only) display cost function surface
    if (displaycf)
        [X,Y] = meshgrid(t_z, r_z);
        figure2 = figure(2);
        surf(X,Y, cf);
        xlabel('tz (mm)');
        ylabel('Rz (deg)');
        zlabel('Cost Function Value');
        view(0, 90);
    end 

    % search min
    minValue = min(cf(:));
    [costfunctions_min(trial, 1), costfunctions_min(trial, 2)] = find(cf == minValue);
    
    % (for debugging only) break the loop
    if (or(displaybone, displaycf))
        break;
    end

end

middle = ceil(loop_z/2);
fprintf('rz_idx: %d\t\t\trt_idx: %d\n', costfunctions_min(1,1)-middle, costfunctions_min(1,2)-middle);
fprintf('rz_est: %.2f deg\trt_est: %.2f mm\n', r_z(costfunctions_min(:, 1)), t_z(costfunctions_min(:, 2))*1000 );
