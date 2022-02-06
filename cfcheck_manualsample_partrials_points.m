clc; clear; close all;

addpath(genpath('..\pointcloudregistration_evaluations'));
addpath(genpath('..\gmmreg\MATLAB'));

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

%% Simulation

noises            = [1 2 3];
pointcounts       = [15 20 25 30];
num_trials        = 150;

costfunction_name  = "gmm";
costfunction_scale = 40;
costfunctions_min  = ones(num_trials, 2, length(noises), length(pointcounts));

for pointcount=1:length(pointcounts)

    current_pointcount = pointcounts(pointcount);
    % Read the simulated a-mode measurement point cloud, which is a subset of Ŭ.
    % These a-mode simulated measurement is manually selected from the bone model.
    selectedpoint_str = sprintf('data/bone/amode_tibia_%d.mat', current_pointcount);
    load(selectedpoint_str);
    U = [ vertcat(amode_prereg1.Position); ...
          vertcat(amode_prereg2.Position); ...
          vertcat(amode_prereg3.Position); ...
          vertcat(amode_mid.Position)]';

    for noise=1:length(noises)

        current_noise = noises(noise);

        for trial=1:num_trials
                
            fprintf('pointcount:%d, noise:%d, trials: %d\n', current_pointcount, current_noise, trial);

            % add isotropic zero-mean gaussian noise to U, simulating noise measurement
            random_noise  = -current_noise/ptCloud_scale + (current_noise/ptCloud_scale + current_noise/ptCloud_scale)*rand(3,size(U, 2));
            model_ptCloud = (U + random_noise)';
            
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
                    
                    % GMM L2 Distance
                    scale = costfunction_scale * 1e-4;
                    [f,~] =  GaussTransform(double(model_ptCloud), double(scene_ptCloud), scale);
                    cf_temp(current_t) = -f;
                    
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

        filename = sprintf('results\\tibia_%s_scale%d_%d.mat', costfunction_name, costfunction_scale, num_trials);
        save(filename, 'costfunctions_min', 'r_z', 't_z');

    % end noises
    end

% end pointcounts
end
