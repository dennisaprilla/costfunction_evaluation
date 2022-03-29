%% Summary
% This code is intended to perform a costfunction examination simulation
% based on curvature similarity function. The idea is, since b-mode produce
% a surface bone curvature (which is basically a 3d curve line), maybe we
% can construct a function which calculates the similarity between the
% expected curvature on the bone model surface and the b-mode measurment

clc; clear; close all;

path_pointcloudregistration = 'D:\Documents\BELANDA\PhD Thesis\Code\MATLAB\pointcloudregistration_evaluations';
path_boneUSsimple           = 'D:\Documents\BELANDA\PhD Thesis\Code\MATLAB\boneUSsimple';

path_bone     = strcat(path_pointcloudregistration, filesep, 'data', filesep, 'bone');
path_bmode    = strcat(path_boneUSsimple, filesep, 'outputs', filesep, 'usmeasurement_b');
path_function1 = strcat(path_boneUSsimple, filesep, 'functions', filesep, 'experimental');
path_function2 = strcat(path_boneUSsimple, filesep, 'functions', filesep, 'geometry');
path_output   = 'results';

addpath(path_bone);
addpath(path_bmode);
addpath(path_function1);
addpath(path_function2);

displaybone = true;
displaycf   = true;
use_boneportion    = true;

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
filename_bmodedata = sprintf('usdata_b_0b');
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

% if use_boneportion is specified, we will use only the portion of the bone
% instead of the whole bone. Portion is obtained from simulation toolbox 
% (refer to file usmeasurement_b.m)
if(use_boneportion)
    % get the bone portion
    U_breve_part = get_boneportion(bmode_simulation.portion, U_breve')';
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

%% Cost Function

scene_ptCloud = U_breve_part';
model_ptCloud = Ub_pointcloud;

[nearest_idx, nearest_dist] = knnsearch(scene_ptCloud, model_ptCloud);















