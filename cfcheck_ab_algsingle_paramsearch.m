clc; clear; close all;

path_pointcloudregistration = '..\pointcloudregistration_evaluations';
path_boneUSsimple           = '..\boneUSsimple';
path_gmmreg                 = '..\gmmreg\MATLAB\GaussTransform';

path_bone      = strcat(path_pointcloudregistration, filesep, 'data', filesep, 'bone');
path_amode     = strcat(path_bone, filesep, filesep, 'amode_accessible_sim2');
path_bmode     = strcat(path_boneUSsimple, filesep, 'outputs', filesep, 'usmeasurement_b');
path_function1 = strcat(path_boneUSsimple, filesep, 'functions', filesep, 'experimental');
path_function2 = strcat(path_boneUSsimple, filesep, 'functions', filesep, 'geometry');
path_output    = 'results';

addpath(path_bone);
addpath(path_amode);
addpath(path_bmode);
addpath(path_gmmreg);
addpath(path_function1);
addpath(path_function2);

displaybone = false;
displaycf   = false;

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
           'Tag', 'plot_bone');
    grid on; axis equal; hold on;
    xlabel('X'); ylabel('Y'); zlabel('Z');
end

%% Prepare the A-mode data
filename_amodedata = sprintf('amode_tibia_15');
filepath_amodedata = sprintf('%s%s%s.mat', path_amode, filesep, filename_amodedata);
load(filepath_amodedata);
Ua_pointcloud = vertcat(amode_all.Position)';

% (for debugging only) show figure for sanity check
if(displaybone)
    grid on; axis equal; hold on;
    plot3( axes1, ...
           Ua_pointcloud(1,:), ...
           Ua_pointcloud(2,:), ...
           Ua_pointcloud(3,:), ...
           'or', 'MarkerFaceColor', 'r', ...
           'Tag', 'plot_bone');
end

%% Prepare the B-mode data
filename_bmodedata = sprintf('usdata_b_2b');
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
           'Tag', 'plot_bone');
end

%% Simulate Search Space

% obtain all combination of z rotation and translation
range = 10;
step  = 0.25;
r_z   = (-range:step:range);
t_z   = (-range/ptCloud_scale:step/ptCloud_scale:range/ptCloud_scale);
% change z rotation to rotation matrix
rs = [ r_z', zeros(length(r_z), 2) ];
Rs = eul2rotm(deg2rad(rs), 'ZYX');
% change z translation to translation vector
ts = [ zeros(2, length(t_z)); t_z];

%% Simulation Setup

% setup the noise constants
noise_level     = 2;
noise_tconst    = 1;
noise_skewconst = 0.025;
noise_Rconst    = 1.25;
% setup the noise
noise_a        = noise_level;
noise_bin_t    = noise_level*noise_tconst;
noise_bin_s    = noise_level*noise_skewconst;
noise_bex_R    = noise_level*noise_Rconst;
noise_bex_t    = noise_level*noise_tconst;

% setup the simulation trials config
costfunction_name       = "gmm";
costfunction_scaleconst = 1e-4;
costfunction_scales_a   = [10 20 30 40 50];
costfunction_scales_b   = [10 20 30 40 50];
costfunction_alphaconst = size(Ua_pointcloud, 1)/size(Ub_pointcloud, 1);
costfunction_alphas     = [0.5 1.0 1.5 2.0];
use_boneportion         = true;
num_trials              = 100;

% if use_boneportion is specified, we will use only the portion of the bone
% instead of the whole bone. Portion is obtained from simulation toolbox 
% (refer to file usmeasurement_b.m)
if(use_boneportion)
    % get the bone portion
    % U_breve_part = get_boneportion(bmode_simulation.portion, U_breve')';
    % U_breve_part = get_boneportion([0.125 0.225], U_breve')';
    U_breve_part = get_boneportion([0.25 0.35], U_breve')';
    % display it
    if(displaybone)
        plot3( axes1, ...
               U_breve_part(1,:), ...
               U_breve_part(2,:), ...
               U_breve_part(3,:), ...
               '.', 'Color', [0.9290 0.6940 0.1250], ...
               'MarkerSize', 0.1, ...
               'Tag', 'plot_bone_full');
    end
end

% naming the filename for result
filepath = 'results\abmode_sim3';

%% Simulation Start

% loop for scale_a
for costfunction_scale_a = 1:length(costfunction_scales_a)
current_scale_a = costfunction_scales_a(costfunction_scale_a);
scale_a = current_scale_a * costfunction_scaleconst;

% loop for scale_b
for costfunction_scale_b = 1:length(costfunction_scales_b)
current_scale_b = costfunction_scales_b(costfunction_scale_b);
scale_b = current_scale_b * costfunction_scaleconst;

% loop for alpha
for costfunction_alpha = 1:length(costfunction_alphas) 
current_alpha = costfunction_alphas(costfunction_alpha);
alpha         = current_alpha * costfunction_alphaconst;

% variable that will contains every global minimum of costfunction
costfunctions_min  = ones(num_trials, 2);
% display the current loop
fprintf('scale_a: %d, scale_b: %d, alpha: %.2f -----------------\n', current_scale_a, current_scale_b, current_alpha);


% for each config, do trials until the number of num_trials
for trial=1:num_trials

    fprintf('trials: %d\n', trial);
    
    % add isotropic zero-mean gaussian noise to U, simulating noise measuremen
    random_noise = -noise_a/ptCloud_scale + (noise_a/ptCloud_scale + noise_a/ptCloud_scale)*rand(3,size(Ua_pointcloud, 2));
    Ua_noised = (Ua_pointcloud + random_noise)';

    % put internal and external noise to b-mode
    Ub_noised = Ub_pointcloud;
    % parameter for internal noise
    t2d_noise = [noise_bin_t/ptCloud_scale 0.5*(noise_bin_t/ptCloud_scale)];
    % t2d_noise = [3/ptCloud_scale 0];
    s_noise   = noise_bin_s;
    [Ub_noised, ~] = bmode_addnoise_internal(Ub_plane, Ub_noised, t2d_noise, s_noise);
    % parameter for external noise
    R_noise   = noise_bex_R;
    t3d_noise = noise_bex_t/ptCloud_scale;
    [Ub_noised, ~] = bmode_addnoise_external(Ub_noised, R_noise, t3d_noise);

    % show figure for sanity check
    if (displaybone)
        grid on; axis equal; hold on;
        plot3( axes1, ...
               Ua_noised(:,1), ...
               Ua_noised(:,2), ...
               Ua_noised(:,3), ...
               'oy', ...
               'Tag', 'plot_bone_full');
        plot3( axes1, ...
               Ub_noised(:,1), ...
               Ub_noised(:,2), ...
               Ub_noised(:,3), ...
               'oy', ...
               'Tag', 'plot_bone_full');
    end

    % prepare variable to contains all costfunction value
    cf  = zeros(length(r_z), length(t_z));

    % prepare a variable which needed for parfor
    loop_z = length(r_z);
    loop_t = length(t_z);

    tic;
    parfor current_z = 1:loop_z

        % temproary variable for parfor
        cf_temp  = zeros(1, loop_t);

        for current_t = 1:loop_t

            % transform Ŭ with respected transformation
            scene_ptCloud_a = ( Rs(:,:,current_z) * U_breve + ts(:,current_t) )';
            scene_ptCloud_b = ( Rs(:,:,current_z) * U_breve_part + ts(:,current_t) )';

            % GMM L2 Distance
            if(strcmp(costfunction_name, "gmm"))
                % amode
                [f_a,~] =  GaussTransform(double(Ua_noised), double(scene_ptCloud_a), scale_a);
                % bmode
                [f_b,~] =  GaussTransform(double(Ub_noised), double(scene_ptCloud_b), scale_b);
                % fusion
                cf_temp(current_t) = -(f_a + (alpha * f_b));

            % RMSE
            elseif (strcmp(costfunction_name, "rmse"))
                % amode
                [~, nearest_dist] = knnsearch(scene_ptCloud_a, Ua_noised);
                f_a = mean(nearest_dist);
                % bmode
                [~, nearest_dist] = knnsearch(scene_ptCloud_b, Ub_noised);
                f_b = mean(nearest_dist);
                % fusion
                cf_temp(current_t) = f_a + (alpha * f_b);
                
            % GMM and RMSE
            elseif (strcmp(costfunction_name, "gmmrmse"))
                % amode
                [f_a,~] =  GaussTransform(double(Ua_noised), double(scene_ptCloud_a), scale_a);
                % bmode
                [~, nearest_dist] = knnsearch(scene_ptCloud_b, Ub_noised);
                f_b = mean(nearest_dist);
                % fusion
                cf_temp(current_t) = -f_a + (alpha * f_b);
                
            % RMSE and GMM
            elseif (strcmp(costfunction_name, "rmsegmm"))
                % amode
                [~, nearest_dist] = knnsearch(scene_ptCloud_a, Ua_noised);
                f_a = mean(nearest_dist);
                % bmode
                [f_b,~] =  GaussTransform(double(Ub_noised), double(scene_ptCloud_b), scale_b);
                % bmode
                cf_temp(current_t) = f_a + (alpha * -f_b);
                
            % end costfunction ifs
            end
            
        % end loop tz
        end

        % store the temp to actual variable contains all costfunction value
        cf(current_z, : )  = cf_temp;
        
	% end loop rz
    end
    toc;
    
    % display cost function surface
    if (displaycf)
        [X,Y] = meshgrid(t_z, r_z);
        figure2 = figure(2);
        surf(X,Y, cf);
        xlabel('tz (mm)');
        ylabel('Rz (deg)');
        zlabel('Cost function');
        view(-90, 90);
    end

    % look for the min
    minValue = min(cf(:));
    [costfunctions_min(trial, 1), costfunctions_min(trial, 2)] = find(cf == minValue);

% end trials    
end

% save the details of the simulation
simulation_config.costfunction_name = costfunction_name;
simulation_config.scaleconst        = costfunction_scaleconst;
simulation_config.scale_a           = current_scale_a;
simulation_config.scale_b           = current_scale_b;
simulation_config.alphaconst        = costfunction_alphaconst;
simulation_config.alpha             = current_alpha;
simulation_config.trials            = trial;
simulation_config.use_boneportion   = use_boneportion;
simulation_config.noise.level       = noise_level;
simulation_config.noise.tconst      = noise_tconst;
simulation_config.noise.skewconst   = noise_skewconst;
simulation_config.noise.Rconst      = noise_Rconst;

% save the necessary files
filename = sprintf('abmodeparamsearch_%d_%d_%d.mat', current_scale_a, current_scale_b, current_alpha*10);
fullpath = strcat(filepath, filesep, filename);
save(fullpath, 'costfunctions_min', 'simulation_config', 'r_z', 't_z');

% end alpha
end
% end scale_b
end
% end scale_a
end
















