clear; close all; clc;
addpath('..\functions\display\subaxis');

% load the data
% filenames={'tibia_gmm_scale20', 'tibia_gmm_scale30', 'tibia_gmm_scale40', 'tibia_rmse'};
filenames={'tibia_gmm_scale40'};

data = [];
for file=1:length(filenames)
    
    current_filename = filenames{file};
    load(strcat(current_filename, '.mat'));

    % setup the noises (this is supposed to be given from the mat file, but i
    % forgot to save it.
    noises            = [1 2 3];
    pointcounts       = [15 20 25 30];

    % calculate the magnitude (maybe it will be used, or not)
    middle = ceil(length(r_z)/2);
    costfunctions_min_normalized = costfunctions_min - middle;
    costfunctions_min_magnitude  = sqrt(sum((costfunctions_min_normalized.^2),2));

    % costfunctions_min contains index of the search-space matrix, let's
    % convert it to real rz and tz value
    tz_scale   = 1000;
    rz_tz_est  = cat(2, r_z(costfunctions_min(:,1,:,:)), t_z(costfunctions_min(:,2,:,:))*tz_scale );

    rz_tz_mean    = mean(rz_tz_est, 1);
    rz_tz_meanabs = mean(abs(rz_tz_est), 1);
    rz_tz_cov     = {};
    
    n_dim       = size(costfunctions_min, 2);
    n_noises    = size(costfunctions_min, 3);
    n_numpoints = size(costfunctions_min, 4);
    for current_noise=1:n_noises
        for current_numpoints=1:n_numpoints
            rz_tz_cov{current_noise, current_numpoints} = diag(cov(rz_tz_est(:,:,current_noise, current_numpoints)))';
        end
    end
    
    %{
    % run this block of lines to understand what am i doing
    A = zeros(1, 2, 3, 4);
    A(:,:,:,1) = [  1  2;  3  4;  5  6 ]';
    A(:,:,:,2) = [  7  8;  9 10; 11 12 ]';
    A(:,:,:,3) = [ 13 14; 15 16; 17 18 ]';
    A(:,:,:,4) = [ 19 20; 21 22; 23 24 ]';
    A
    A_squeezed = squeeze(A)
    A_reshaped = reshape(A_squeezed, 2, 12)'
    %}
    allmean    = reshape(squeeze(rz_tz_mean), n_dim, n_noises*n_numpoints)';
    allmeanabs = reshape(squeeze(rz_tz_meanabs), n_dim, n_noises*n_numpoints)';
    
    %{
    % run this block of lines to understand what am i doing
    A = [ 1 2 3; 4 5 6; 7 8 9; 10 11 12 ]'
    A_reshaped = reshape(A, 12, 1)
    %}
    alldiagcov = cell2mat(reshape(rz_tz_cov, 12, 1));
    
    % put everything to data variable
    data = [data, allmeanabs, allmean, alldiagcov];

end

data = round(data, 3);

