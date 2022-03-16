clc; clear; close all;

path_pointcloudregistration = 'D:\Documents\BELANDA\PhD Thesis\Code\MATLAB\pointcloudregistration_evaluations';
path_boneUSsimple           = 'D:\Documents\BELANDA\PhD Thesis\Code\MATLAB\boneUSsimple';
path_gmmreg                 = 'D:\Documents\BELANDA\PhD Thesis\Code\cpp\gmmreg\MATLAB\GaussTransform';

path_bone     = strcat(path_pointcloudregistration, filesep, 'data', filesep, 'bone');
path_bmode    = strcat(path_boneUSsimple, filesep, 'outputs', filesep, 'usmeasurement_b');
path_function1 = strcat(path_boneUSsimple, filesep, 'functions', filesep, 'experimental');
path_function2 = strcat(path_boneUSsimple, filesep, 'functions', filesep, 'geometry');
path_output   = 'results';

addpath(path_bone);
addpath(path_bmode);
addpath(path_gmmreg);
addpath(path_function1);
addpath(path_function2);

displaybone = true;
displaycf   = true;

clear path_pointcloudregistration path_boneUSsimple path_gmmreg;

%% Prepare Bone Point Cloud

% read the point cloud (bone) from STL/PLY file
ptCloud          = stlread('data/bone/CT_Tibia_R.stl');
ptCloud_scale    = 1000;
ptCloud_Npoints  = size(ptCloud.Points,1);
ptCloud_centroid = mean(ptCloud.Points, 1);
% prepare Ŭ, the noiseless, complete, moving dataset
U_breve          = (ptCloud.Points - ptCloud_centroid)';

% show figure for sanity check
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
    grid on; axis equal; hold on;
    xlabel('X'); ylabel('Y'); zlabel('Z');
end


%% Prepare the B-mode data
filename_bmodedata = sprintf('usdata_b_1a');
filepath_bmodedata = sprintf('%s%s%s.mat', path_bmode, filesep, filename_bmodedata);
load(filepath_bmodedata);
Ub_pointcloud = bmode_simulation.pointcloud;
Ub_plane      = bmode_simulation.plane;

% show figure for sanity check
if(displaybone)
    plot3( axes1, ...
           Ub_pointcloud(:,1), ...
           Ub_pointcloud(:,2), ...
           Ub_pointcloud(:,3), ...
           'or', 'MarkerFaceColor', 'r', ...
           'Tag', 'plot_bone_full');
end

%% Simulate Search Space

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

%% Simulation Setup

% setup the noise constants
noise_level     = 1;
noise_skewconst = 0.025;
noise_Rconst    = 1.50;
% setup the noise
noise_bin_t    = noise_level;
noise_bin_s    = noise_level*noise_skewconst;
noise_bex_R    = noise_level*noise_Rconst;
noise_bex_t    = noise_level;

% setup the simulation trials config
num_trials         = 1;
costfunction_name  = "gmm";
costfunction_scale = 10;
use_boneportion    = true;

% if use_boneportion is specified, we will use only the portion of the bone
% instead of the whole bone. Portion is obtained from simulation toolbox 
% (refer to file usmeasurement_b.m)
if(use_boneportion)
    % get the bone portion
    U_breve = get_boneportion(bmode_simulation.portion, U_breve')';
    % display it
    if(displaybone)
        plot3( axes1, ...
               U_breve(1,:), ...
               U_breve(2,:), ...
               U_breve(3,:), ...
               '.', 'Color', [0.9290 0.6940 0.1250], ...
               'MarkerSize', 0.1, ...
               'Tag', 'plot_bone_full');
    end
end

% variable that will contains every global minimum of costfunction
costfunctions_min  = ones(num_trials, 2);
% naming the filename for result
filename_simresult = sprintf('tibia_ab_%s_%d', costfunction_name, num_trials);
    
% for each config, do trials until the number of num_trials
for trial=1:num_trials

    fprintf('trials: %d\n', trial);

    % put internal and external noise to b-mode
    Ub_noised = Ub_pointcloud;
    % parameter for internal noise
    t2d_noise = [noise_bin_t/ptCloud_scale 0.5*(noise_bin_t/ptCloud_scale)];
    s_noise   = noise_bin_s;
    [Ub_noised, ~] = bmode_addnoise_internal(Ub_plane, Ub_noised, t2d_noise, s_noise);
    % parameter for external noise
    R_noise   = noise_bex_R;
    t3d_noise = noise_bex_t/ptCloud_scale;
    [Ub_noised, ~] = bmode_addnoise_external(Ub_noised, R_noise, t3d_noise);

    % rename variable
    model_ptCloud = Ub_noised;

    % show figure for sanity check
    if (displaybone)
        grid on; axis equal; hold on;
        plot3( axes1, ...
               model_ptCloud(:,1), ...
               model_ptCloud(:,2), ...
               model_ptCloud(:,3), ...
               'oy', ...
               'Tag', 'plot_bone_full');
    end

    % prepare variable to contains all costfunction value
    cf  = zeros(length(r_z), length(t_z));

    % prepare a variable which needed for par for
    loop_z = length(r_z);
    loop_t = length(t_z);

    tic;
    parfor current_z = 1:loop_z

        % temproary variable for parfor
        cf_temp  = zeros(1, loop_t);

        for current_t = 1:loop_t

            % transform Ŭ with respected transformation
            U_breve_prime = Rs(:,:,current_z) * U_breve + ts(:,current_t);
            scene_ptCloud = U_breve_prime';

            if(strcmp(costfunction_name, "gmm"))
                % GMM L2 Distance
                scale = costfunction_scale * 1e-4;
                [f,~] =  GaussTransform(double(model_ptCloud), double(scene_ptCloud), scale);
                cf_temp(current_t) = -f;

            elseif (strcmp(costfunction_name, "rmse"))
                % RMSE
                [nearest_idx, nearest_dist] = knnsearch(scene_ptCloud, model_ptCloud);
                cf_temp(current_t) = mean(nearest_dist);
            end

        end

        % store the temp to actual variable contains all costfunction value
        cf(current_z, : )  = cf_temp;
    end
    toc;
    
    % display cost function surface
    if (displaycf)
        [X,Y] = meshgrid(t_z, r_z);
        figure2 = figure(2);
        surf(X,Y, cf);
        xlabel('tz (mm)');
        ylabel('Rz (deg)');
        zlabel('GMM L2 distance');
        view(-90, 90);
    end

    % look for the min
    minValue = min(cf(:));
    [costfunctions_min(trial, 1), costfunctions_min(trial, 2)] = find(cf == minValue);

% end trials    
end
    
middle = ceil(loop_z/2);
fprintf('rz_idx: %d\t\t\trt_idx: %d\n', costfunctions_min(1,1)-middle, costfunctions_min(1,2)-middle);
fprintf('rz_est: %.2f deg\trt_est: %.2f mm\n', r_z(costfunctions_min(:, 1)), t_z(costfunctions_min(:, 2))*1000 );
















