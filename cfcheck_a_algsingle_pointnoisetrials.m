clc; clear; close all;

addpath(genpath('..\pointcloudregistration_evaluations'));
addpath(genpath('..\gmmreg\MATLAB'));

displaybone = true;
displaycf   = true;

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

%% Simulation

% setup the simulation configuration
noises             = [1 2 3];
pointcounts        = [10 15 20 25 30];
num_trials         = 500;
costfunction_name  = "rmse";
costfunction_scale = 20;

% save the configuration to a structure
trialsdesc.noises              = noises;
trialsdesc.pointconfigs        = pointconfigs;
trialsdesc.num_trials          = num_trials;
trialsdesc.costfunction_name   = costfunction_name;
if (strcmp(costfunction_name, "gmm"))
    trialsdesc.costfunction_scale  = costfunction_scale;
    % naming the filename for result
    filename_simresult = sprintf('tibia_%s_scale%d_%d', costfunction_name, costfunction_scale, num_trials);
else
    trialsdesc.costfunction_scale  = NaN;
    % naming the filename for result
    filename_simresult = sprintf('tibia_%s_%d', costfunction_name, num_trials);
end

% variable that will contains every global minimum of costfunction
costfunctions_min  = ones(num_trials, 2, length(noises), length(pointcounts));

% loop over all of the pointcount configuration
for pointconfig=1:length(pointconfigs)

    current_pointconfig = pointconfigs(pointconfig);
    
    % Read the simulated a-mode measurement point cloud, which is a subset of Ŭ.
    % These a-mode simulated measurement is manually selected from the bone model.
    filename_amodedata = sprintf('amode_tibia_%d', current_pointconfig);
    filepath_amodedata = sprintf('data/bone/amode_accessible_sim2/%s.mat', filename_amodedata);
    load(filepath_amodedata);
    U = [ vertcat(amode_prereg1.Position); ...
          vertcat(amode_prereg2.Position); ...
          vertcat(amode_prereg3.Position); ...
          vertcat(amode_mid.Position)]';
    %}
    %
    filename_amodedata = sprintf('amode_tibia_%d', current_pointcount);
    filepath_amodedata = sprintf('data/bone/amode_accessible_sim2/%s.mat', filename_amodedata);
    load(filepath_amodedata);
    U = vertcat(amode_all.Position)';
    %
      
    % loop over all of the noise configuration
    for noise=1:length(noises)

        current_noise = noises(noise);
        
        % for each config, do trials until the number of num_trials
        for trial=1:num_trials
                
            fprintf('pointcount:%d, noise:%d, trials: %d\n', current_pointconfig, current_noise, trial);

            % add isotropic zero-mean gaussian noise to U, simulating noise measurement
            random_noise  = -current_noise/ptCloud_scale + (current_noise/ptCloud_scale + current_noise/ptCloud_scale)*rand(3,size(U, 2));
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
                zlabel('GMM L2 distance');
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
        save(sprintf('results\\accessible_sim2\\%s.mat', filename_simresult), 'costfunctions_min', 'r_z', 't_z', 'trialsdesc');

    % end noises
    end
    
    % (for debugging only) break the nested loop
    if (or(displaybone, displaycf))
        break;
    end

% end pointcounts
end
