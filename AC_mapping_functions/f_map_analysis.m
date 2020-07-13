function f_map_analysis(data, params, ops)

disp('Analysis...');

vec_frame_data = data.vec_frame_data;
trial_types = data.trial_types;
%volt_data_binned = data.volt_data_binned;
clear data;

[d1, d2, num_frames, num_trials] = size(vec_frame_data);

fig1 = figure;
imagesc(squeeze(mean(mean(vec_frame_data,3),4)));
colormap gray;
title(['Average Image ' ops.file_name], 'Interpreter', 'none');
axis tight equal off;
savefig(fig1, [ops.data_dir '\Average Image ' ops.file_name]);

%% normalize pizel by pixel

vec_frame_data = f_map_norm_by_pix(vec_frame_data);
vec_frame_data = single(reshape(vec_frame_data, d1, d2, num_frames, num_trials));

% figure;
% imagesc(squeeze((vec_frame_data(:,:,1,1))));
% colormap gray;
% title(['Average Image norm ' ops.file_name], 'Interpreter', 'none');
% axis tight equal off;


% save something, if needed to check normalization

% size(norm2_vec_data(:,:,:,1:50))
% vid_stack = reshape(norm2_vec_data(:,:,:,1:50), [res1 res2 sta_frames * 50]);
% imwrite((vid_stack(:,:,1)), sprintf('%s_normalized_test.tif',file_name))
% for jj = 2:size(vid_stack,3)
%     imwrite((vid_stack(:,:,jj)), sprintf('%s_normalized_test.tif',file_name), 'WriteMode', 'append');
% end
% 

%%
disp('Plotting...');
figure;
histogram(trial_types,params.num_trial_types);
title(sprintf('stimulus counts, total = %d', length(trial_types)));

%% trial averages
f_map_plot_trial_ave(vec_frame_data, trial_types, params, ops);

%%
f_map_plot_spatial_trial_ave(vec_frame_data, trial_types, params, ops)

%% analyze single trials
if ops.save_trial_analysis
    f_map_single_trial_analysis(vec_frame_data, trial_types, params, ops)
end


%% trial ave vids
if ops.save_trial_ave_vids
    f_map_trial_ave_vid(vec_frame_data, trial_types, params, ops)
end

%% concat single trials
if ops.save_single_trial_vids
    f_map_single_trial_vid(vec_frame_data, trial_types, params, ops)
end

end