clc; clear; close all;

path_pointcloudregistration = 'D:\Documents\BELANDA\PhD Thesis\Code\MATLAB\pointcloudregistration_evaluations';
path_boneUSsimple           = 'D:\Documents\BELANDA\PhD Thesis\Code\MATLAB\boneUSsimple';
path_gmmreg                 = 'D:\Documents\BELANDA\PhD Thesis\Code\cpp\gmmreg\MATLAB';

path_bone     = strcat(path_pointcloudregistration, filesep, 'data', filesep, 'bone');
path_amode    = strcat(path_bone, filesep, 'amode_accessible_sim3');
path_bmode    = strcat(path_boneUSsimple, filesep, 'outputs', filesep, 'usmeasurement_b');
path_function = strcat(path_boneUSsimple, filesep, 'functions');
path_output   = 'results';

addpath(path_bone);
addpath(path_amode);
addpath(path_bmode);
addpath(genpath(path_function));

clear path_pointcloudregistration path_boneUSsimple path_gmmreg;

%% Prepare Bone Point Cloud

% read the point cloud (bone) from STL/PLY file
ptCloud          = stlread('data/bone/CT_Tibia_R.stl');
ptCloud_scale    = 1000;
ptCloud_Npoints  = size(ptCloud.Points,1);
ptCloud_centroid = mean(ptCloud.Points, 1);
% prepare Ŭ, the noiseless, complete, moving dataset
U_breve          = (ptCloud.Points - ptCloud_centroid)';

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

% setup the simulation configuration
noises             = [1 2];
pointcounts        = 15;
num_trials         = 5;
costfunction_name  = "rmse";
costfunction_scale = 40;

% save the configuration to a structure
trialsdesc.noises              = noises;
trialsdesc.pointcounts         = pointcounts;
trialsdesc.num_trials          = num_trials;
trialsdesc.costfunction_name   = costfunction_name;
if (strcmp(costfunction_name, "gmm"))
    trialsdesc.costfunction_scale  = costfunction_scale;
else
    trialsdesc.costfunction_scale  = NaN;
end

% variable that will contains every global minimum of costfunction
costfunctions_min  = ones(num_trials, 2, length(noises), length(pointcounts));
% naming the filename for result
filename_simresult = sprintf('tibia_ab_%s_%d', costfunction_name, num_trials);

%%

% loop over all of the pointcount configuration
for pointcount=1:length(pointcounts)

    current_pointcount = pointcounts(pointcount);
    
    % get a-mode data
    filename_amodedata = sprintf('amode_tibia_%d_a', current_pointcount);
    filepath_amodedata = sprintf('%s%s%s.mat', path_amode, filesep, filename_amodedata);
    load(filepath_amodedata);
    Ua = vertcat(amode_all.Position);
    
    % get b-mode data
    filename_bmodedata = sprintf('usdata_b_02-23-2022_11-33');
    filepath_bmodedata = sprintf('%s%s%s.mat', path_bmode, filesep, filename_bmodedata);
    load(filepath_bmodedata);
    Ub_pointcloud = bmode_simulation.pointcloud;
    Ub_plane      = bmode_simulation.plane;
    
    %{
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
    grid on; axis equal; hold on;
    plot3( axes1, ...
           Ua(1,:), ...
           Ua(2,:), ...
           Ua(3,:), ...
           'or', 'MarkerFaceColor', 'r', ...
           'Tag', 'plot_bone_full');
    plot3( axes1, ...
           Ub(1,:), ...
           Ub(2,:), ...
           Ub(3,:), ...
           'or', 'MarkerFaceColor', 'r', ...
           'Tag', 'plot_bone_full');
    %}
    
    % loop over all of the noise configuration
    for noise=1:length(noises)

        current_noise = noises(noise);
        
        % for each config, do trials until the number of num_trials
        for trial=1:num_trials
            
            fprintf('pointcount:%d, noise:%d, trials: %d\n', current_pointcount, current_noise, trial);
            
            % add isotropic zero-mean gaussian noise to U, simulating noise measurement
            random_noise  = -current_noise/ptCloud_scale + (current_noise/ptCloud_scale + current_noise/ptCloud_scale)*rand(size(Ua, 1),3);
            Ua_noised = (Ua + random_noise);
            
            % put internal and external noise to b-mode
            t2d_noise = current_noise/ptCloud_scale;
            s_noise   = current_noise/2 * 0.1;
            [Ub_noised, ~] = bmode_addnoise_internal(Ub_plane, Ub_pointcloud, t2d_noise, s_noise);
            R_noise   = current_noise;
            t3d_noise = current_noise/ptCloud_scale;
            [Ub_noised, ~] = bmode_addnoise_external(Ub_noised, R_noise, t3d_noise);            

            % gather all of the ultrasound measurement simulation
            model_ptCloud = [Ua_noised; Ub_noised];
            
            %
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
            grid on; axis equal; hold on;
            plot3( axes1, ...
                   model_ptCloud(:,1), ...
                   model_ptCloud(:,2), ...
                   model_ptCloud(:,3), ...
                   'or', 'MarkerFaceColor', 'r', ...
                   'Tag', 'plot_bone_full');
            %
            
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
            
            % look for the min
            minValue = min(cf(:));
            [costfunctions_min(trial, 1, noise, pointcount), ...
             costfunctions_min(trial, 2, noise, pointcount)] = find(cf == minValue);
         
        % end trials    
        end
        
        % i put save here, just in case the pc is overheating
        save(sprintf('%s%s%s.mat', path_output, filesep, filename_simresult), 'costfunctions_min', 'r_z', 't_z', 'trialsdesc');
    
    % end noises    
    end
    
% end pointcounts    
end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
