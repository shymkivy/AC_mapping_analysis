function f_map_trial_ave_vid(vec_frame_data, trial_types, params, ops)

concat_vids = 1;
add_sig_indicator = 1;

[d1, d2, num_frames, ~] = size(vec_frame_data);
freq_amp_lookup = params.freq_amp_lookup;
stim_type_lookup = params.stim_type_lookup;

if ~exist([ops.data_dir '\trial_ave_vids'], 'dir')
    mkdir([ops.data_dir '\trial_ave_vids']);
end


ave_frames = zeros(d1, d2, num_frames, params.num_trial_types);
for n_tr = 1:params.num_trial_types
    ave_frames(:,:,:,n_tr) = mean(vec_frame_data(:,:,:,trial_types == n_tr),4);
end

if concat_vids
    m = 2; n = 5;
    for n_amp = 1:numel(params.stim_params.modulation_amp)
        temp_ave_data = ave_frames(:,:,:,freq_amp_lookup(:, n_amp));
        % concatenate spatially
        temp_ave_data1 = reshape(permute(temp_ave_data, [2 4 1 3]), d2*n,[],d1,num_frames);
        temp_ave_data2 = reshape(permute(temp_ave_data1, [3 2 1 4]), d1*m, d2*n, num_frames);
        
        
        vid_stack = temp_ave_data2;
        vid_stack = vid_stack - mean(vid_stack(:,:,params.baseline_window_frames),3);
        if add_sig_indicator
            sig_trace = ones(5,size(vid_stack,2), size(vid_stack,3))*min(vid_stack(:));
            sig_trace(:,:,params.onset_window_frames) = max(vid_stack(:));
            vid_stack = cat(1, sig_trace, vid_stack);
        end
        

        f_save_tif_stack2_YS(vid_stack, sprintf('%s_allfr_%dV_trial_ave.tif',[ops.data_dir '\trial_ave_vids\' ops.file_name], stim_type_lookup(freq_amp_lookup(1, n_amp),2)));
        %f_play_stack(vid_stack, 30, 1, [num2str(all_trial_types(freq_amp_lookup(1, n_amp),2)) 'V']);
        
    end
else
    for n_amp = 1:numel(params.stim_params.modulation_amp)
        for n_fr = 1:params.stim_params.num_freqs
            vid_stack = squeeze(ave_frames(:,:,:,freq_amp_lookup(n_fr, n_amp)));
            vid_stack = vid_stack - mean(vid_stack(:,:,params.baseline_window_frames),3);

            f_save_tif_stack_YS(vid_stack, sprintf('%s_%dV_%.1fkHz_trial_average.tif',[ops.data_dir '\trial_ave_vids\' ops.file_name], stim_type_lookup(freq_amp_lookup(n_fr, n_amp),2), stim_type_lookup(freq_amp_lookup(n_fr, n_amp),1)/1000));
            %implay(vid_stack);
            %implay(uint8(vid_stack*(2^8)), 20);
        end
    end
end

end