function f_map_single_trial_vid(vec_frame_data, trial_types, params, ops)

[d1, d2, num_frames, ~] = size(vec_frame_data);
freq_amp_lookup = params.freq_amp_lookup;

if ~exist([ops.data_dir '\single_trial_vids'], 'dir')
    mkdir([ops.data_dir '\single_trial_vids']);
end


for n_amp = 1:numel(params.stim_params.modulation_amp)
    for n_fr = 1:params.stim_params.num_freqs  
        trial_vec_data = vec_frame_data(:,:,:,trial_types == freq_amp_lookup(n_fr, n_amp));
        vid_stack = f_make_single_trial_tiles(trial_vec_data, params.onset_window_frames); 
        f_save_tif_stack2_YS(vid_stack, sprintf('%s_%dV_%.1fkHz_single_trials.tif',[ops.data_dir '\single_trial_vids\' ops.file_name], params.stim_type_lookup(freq_amp_lookup(n_fr, n_amp),2), params.stim_type_lookup(freq_amp_lookup(n_fr, n_amp),1)/1000));
        %f_play_stack(vid_stack, 30, 1, [num2str(params.stim_type_lookup(freq_amp_lookup(n_fr, n_amp),1)/1000) ' kHz']);
    end
end


end