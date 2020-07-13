function f_map_single_trial_vid(vec_frame_data, trial_types, params, ops)

add_sig_indicator = 1;

[d1, d2, num_frames, ~] = size(vec_frame_data);
freq_amp_lookup = params.freq_amp_lookup;

if ~exist([ops.data_dir '\single_trial_vids'], 'dir')
    mkdir([ops.data_dir '\single_trial_vids']);
end


for n_amp = 1:numel(params.stim_params.modulation_amp)
    for n_fr = 1:params.stim_params.num_freqs
        
        trial_vec_data = vec_frame_data(:,:,:,trial_types == freq_amp_lookup(n_fr, n_amp));
        num_trials = sum(trial_types == freq_amp_lookup(n_fr, n_amp));
        n = 10;
        m = ceil(num_trials/n);
        
        % concatenate spatially
        temp_fill = ones(d1, d2, num_frames, n*m-num_trials)*min(trial_vec_data(:));
        trial_vec_data = cat(4, trial_vec_data, temp_fill);
        trial_vec_data1 = reshape(permute(trial_vec_data, [2 4 1 3]), d2*n,[],d1,num_frames);
        trial_vec_data2 = reshape(permute(trial_vec_data1, [3 2 1 4]), d1*m, d2*n, num_frames);

        vid_stack = trial_vec_data2;
        %vid_stack = vid_stack - mean(vid_stack(:,:,params.baseline_window_frames),3);
        %vid_stack = (vid_stack - min(vid_stack(:)))/(max(vid_stack(:)) - min(vid_stack(:)));
        if add_sig_indicator
            sig_trace = ones(5,size(vid_stack,2), size(vid_stack,3))*min(vid_stack(:));
            sig_trace(:,:,params.onset_window_frames) = max(vid_stack(:));
            vid_stack = cat(1, sig_trace, vid_stack);
        end
        
        f_save_tif_stack_YS(vid_stack, sprintf('%s_%dV_%.1fkHz_single_trials.tif',[ops.data_dir '\single_trial_vids\' ops.file_name], params.stim_type_lookup(freq_amp_lookup(n_fr, n_amp),2), params.stim_type_lookup(freq_amp_lookup(n_fr, n_amp),1)/1000));
        %f_play_stack(vid_stack, 30, 1, [num2str(params.stim_type_lookup(freq_amp_lookup(n_fr, n_amp),1)/1000) ' kHz']);
    end
end


end