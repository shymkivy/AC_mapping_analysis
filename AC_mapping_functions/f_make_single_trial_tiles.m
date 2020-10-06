function vid_stack = f_make_single_trial_tiles(trial_vec_data, stim_frames)

[d1,d2,num_frames,num_trials] = size(trial_vec_data);

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
if exist('stim_frames', 'var')
    sig_trace = ones(5,size(vid_stack,2), size(vid_stack,3))*min(vid_stack(:));
    sig_trace(:,:,stim_frames) = max(vid_stack(:));
    vid_stack = cat(1, sig_trace, vid_stack);
end


end