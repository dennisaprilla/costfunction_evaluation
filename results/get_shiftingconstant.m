clear; close all; clc;

filepath  = 'amode_simulations\accessible_sim2';
filenames = {'tibia_gmm_scale40'};

for file=1:length(filenames)
    
    current_filename = filenames{file};
    load(strcat(filepath, filesep, current_filename, '.mat'));
    
    % costfunctions_min contains index of the search-space matrix, let's
    % convert it to real rz and tz value
    rz_tz_est  = cat(2, r_z(costfunctions_min(:,1,:,:)), t_z(costfunctions_min(:,2,:,:)) );
    
    rz_tz_mean    = mean(rz_tz_est, 1);
    sc            = mean(squeeze(rz_tz_mean),2);
    
    keyset = num2cell( trialsdesc.pointcounts);
    valset = {};
    for pointcount=1:length(trialsdesc.pointcounts)
        valset{pointcount} = sc(:,:,pointcount);
    end
    shiftingconstant = containers.Map(keyset, valset);
    
    save( strcat(filepath, filesep, current_filename, '_shiftingconstant.mat'), ...
         'shiftingconstant');
    
end

