clc; clear; close all;

path_pointcloudregistration = '..\pointcloudregistration_evaluations';
path_boneUSsimple           = '..\boneUSsimple';
path_gmmreg                 = '..\gmmreg\MATLAB\GaussTransform';

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

% set both of these to false if you are not using debug mode
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
noises             = [1 2 3];
noise_skewconst    = 0.025;
noise_Rconst       = 1.5;
pointconfigs       = {'usdata_b_1a', 'usdata_b_1b', 'usdata_b_2a', 'usdata_b_2b'};
num_trials         = 500;
costfunction_name  = "gmm";
costfunction_scale = 20;

% save the configuration to a structure
trialsdesc.noises              = noises;
trialsdesc.pointconfigs        = pointconfigs;
trialsdesc.num_trials          = num_trials;
trialsdesc.costfunction_name   = costfunction_name;
if (strcmp(costfunction_name, "gmm"))
    trialsdesc.costfunction_scale  = costfunction_scale;
else
    trialsdesc.costfunction_scale  = NaN;
end

% variable that will contains every global minimum of costfunction
costfunctions_min  = ones(num_trials, 2, length(noises), length(pointconfigs));
% naming the filename for result
filename_simresult = sprintf('tibiaBmode_%s%d_%d', trialsdesc.costfunction_name, trialsdesc.costfunction_scale, trialsdesc.num_trials);

%%

% loop over all of the pointconfig configuration
for pointconfig=1:length(pointconfigs)
    
    % get the current pointconfig
    current_pointconfig = pointconfigs{pointconfig};
    
    % get b-mode data
    filename_bmodedata = sprintf(current_pointconfig);
    filepath_bmodedata = sprintf('%s%s%s.mat', path_bmode, filesep, filename_bmodedata);
    load(filepath_bmodedata);
    Ub_pointcloud = bmode_simulation.pointcloud;
    Ub_plane      = bmode_simulation.plane;
    
    % (for debugging only) show figure for sanity check
    if(displaybone)
        grid on; axis equal; hold on;
        plot3( axes1, ...
               Ub_pointcloud(:,1), ...
               Ub_pointcloud(:,2), ...
               Ub_pointcloud(:,3), ...
               'or', 'MarkerFaceColor', 'r', ...
               'Tag', 'plot_bone_full');
    end
    
    % loop over all of the noise configuration
    for noise=1:length(noises)
        
        % get the current noise
        current_noise = noises(noise);
        
        % setup the noise
        noise_bin_t    = current_noise;
        noise_bin_s    = current_noise*noise_skewconst;
        noise_bex_R    = current_noise*noise_Rconst;
        noise_bex_t    = current_noise;
        
        % for each config, do trials until the number of num_trials
        for trial=1:num_trials
            
            fprintf('pointconfig: %s, noise: %d, trials: %d\n', current_pointconfig, current_noise, trial);
            
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

            % gather all of the ultrasound measurement simulation
            model_ptCloud = Ub_noised;
            
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
            
            % look for the min
            minValue = min(cf(:));
            [costfunctions_min(trial, 1, noise, pointconfig), ...
             costfunctions_min(trial, 2, noise, pointconfig)] = find(cf == minValue);
         
            % (for debugging only) break the nested loop
            if (or(displaybone, displaycf))
                break;
            end
         
        % end trials    
        end
        
        % i put save here, just in case the pc is overheating
        save(sprintf('%s%s%s.mat', path_output, filesep, filename_simresult), 'costfunctions_min', 'r_z', 't_z', 'trialsdesc');
        
        % (for debugging only) break the nested loop
        if (or(displaybone, displaycf))
            break;
        end
    
    % end noises    
    end
    
    % (for debugging only) break the nested loop
    if (or(displaybone, displaycf))
        break;
    end
    
% end pointconfigs    
end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
