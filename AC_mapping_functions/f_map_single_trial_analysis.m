function f_map_single_trial_analysis(vec_frame_data, trial_types, params, ops)

[d1, d2, num_frames, ~] = size(vec_frame_data);
freq_amp_lookup = params.freq_amp_lookup;


for n_amp = 1:numel(params.stim_params.modulation_amp)
    for n_fr = 2:params.stim_params.num_freqs
        
        temp_vec_data = vec_frame_data(:,:,:,trial_types == freq_amp_lookup(n_fr, n_amp));
        num_trials = size(temp_vec_data,4);
        
        kernel_size = 20;
        kernel = hamming(kernel_size)*hamming(kernel_size)';
        kernel = kernel/sum(kernel(:));

        temp_vec_data_sm = zeros(size(temp_vec_data));
        for n_fr2 = 1:num_frames
            for n_tr = 1:num_trials
                temp_vec_data_sm(:,:,n_fr2, n_tr) = conv2(temp_vec_data(:,:,n_fr2,n_tr),kernel, 'same');
            end
        end
        temp_vec_data_sm = imresize(temp_vec_data_sm, 0.25);
        
        %% save vid
        [d1s, d2s, ~, ~] = size(temp_vec_data_sm);
        n = 10;
        m = ceil(num_trials/n);
        % concatenate spatially
        temp_fill = ones(d1s, d2s, num_frames, n*m-num_trials)*min(temp_vec_data_sm(:));
        temp_vec_data1 = cat(4, temp_vec_data_sm, temp_fill);
        temp_vec_data2 = reshape(permute(temp_vec_data1, [2 4 1 3]), d2s*n,[],d1s,num_frames);
        temp_vec_data3 = reshape(permute(temp_vec_data2, [3 2 1 4]), d1s*m, d2s*n, num_frames);

        vid_stack = temp_vec_data3;
        vid_stack = vid_stack - mean(vid_stack(:,:,params.baseline_window_frames),3);
        %vid_stack = (vid_stack - min(vid_stack(:)))/(max(vid_stack(:)) - min(vid_stack(:)));

        sig_trace = ones(5,size(vid_stack,2), size(vid_stack,3))*min(vid_stack(:));
        sig_trace(:,:,params.onset_window_frames) = max(vid_stack(:));
        vid_stack = cat(1, sig_trace, vid_stack);
        temp_frames = logical(params.onset_window_frames+params.baseline_window_frames);
        %vid_stack = vid_stack(:,:,temp_frames,:);
        
        f_save_tif_stack_YS(vid_stack, sprintf('%s_%dV_%.1fkHz_test_vid.tif',[ops.data_dir '\' ops.file_name], params.stim_type_lookup(freq_amp_lookup(n_fr, n_amp),2), params.stim_type_lookup(freq_amp_lookup(n_fr, n_amp),1)/1000));
        
        
        
        %% dim reduction

        trial_vec_data_2d = reshape(temp_vec_data_sm, d1s*d2s,[]);
        trial_means = mean(trial_vec_data_2d,2);
        
        num_comp = 5;
 
        
%         [icasig, A, W] = fastica(trial_vec_data_2d,'lastEig', num_comp);
%         factors = A;

%         [d_W,d_H] = nnmf(trial_vec_data_2d,num_comp);
%         factors = d_W;

        est_factors = cp_als(tensor(reshape(trial_vec_data_2d-trial_means, d1s*d2s,[], num_trials)),num_comp);
        factors = est_factors.U{1};
        
        for n_comp = 1:num_comp
            figure;
            imagesc(reshape(factors(:,n_comp),d1s,d2s));
            axis equal tight
            title(sprintf('Feq%d Comp %d',n_fr, n_comp));
        end
        


    end
end

end