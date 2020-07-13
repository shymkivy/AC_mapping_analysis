function f_map_plot_spatial_trial_ave(vec_frame_data, trial_types, params, ops)
norm_axis = 0;

[d1, d2, num_frames, num_trials] = size(vec_frame_data);

%% select trials
freq_amp_lookup = params.freq_amp_lookup;
modulation_amp = params.stim_params.modulation_amp;
num_freqs = params.stim_params.num_freqs;
stim_type_lookup = params.stim_type_lookup;

%% define windows

ave_frames = zeros(d1, d2, num_frames, params.num_trial_types);
for n_tr = 1:params.num_trial_types
    ave_frames(:,:,:,n_tr) = mean(vec_frame_data(:,:,:,trial_types == n_tr),4);
end

for n_amp = 1:numel(modulation_amp)
    c_min = 0;
    c_max = 0;
    fig1 = figure;
    sp = cell(num_freqs,1);
    for n_fr=1:num_freqs
        sp{n_fr} = f_subplot_tight(ops.subplot_dimensions(1),ops.subplot_dimensions(2),n_fr);
        tmp_frame = mean(ave_frames(:,:,params.onset_window_frames,freq_amp_lookup(n_fr,n_amp)),3);
        % subtract baseline if you want
        tmp_frame = tmp_frame - mean(ave_frames(:,:,params.baseline_window_frames,freq_amp_lookup(n_fr,n_amp)),3);
        if ops.smooth_frames
            tmp_frame = conv2(tmp_frame,kernel, 'same');
        end
        c_min = min([c_min; tmp_frame(:)]);
        c_max = max([c_max; tmp_frame(:)]);
        imagesc(tmp_frame);
        %colormap gray
        axis equal tight off;
        
        title(['\fontsize{10}' sprintf('diff %d: %.1f kHz', n_fr, stim_type_lookup(n_fr)/1000)]);
    end
    if norm_axis
        for n_fr=1:num_freqs
            caxis(sp{n_fr}, [c_min c_max]);
        end
    end
    suptitle(sprintf('%dV modulation amplitude', modulation_amp(n_amp)));
end
savefig(fig1, [ops.data_dir '\Spatial_ave ' ops.file_name]);

%% combine over amplitudes
if numel(modulation_amp)>1
    ave_frames = zeros(d1, d2, num_frames, num_freqs);
    for n_tr = 1:num_freqs
        ave_frames(:,:,:,n_tr) = mean(vec_frame_data(:,:,:,logical(sum(trial_types == freq_amp_lookup(n_tr,:),2))),4);
    end

    fig1 = figure;
    for n_fr=1:num_freqs
        subplot(ops.subplot_dimensions(1),ops.subplot_dimensions(2),n_fr);
        tmp_frame = mean(ave_frames(:,:,params.onset_window_frames,n_fr),3);
        % subtract baseline if you want
        tmp_frame = tmp_frame - mean(ave_frames(:,:,params.baseline_window_frames,n_fr),3);
        if ops.smooth_frames
            tmp_frame = conv2(tmp_frame,kernel, 'same');
        end
        imagesc(tmp_frame);
        %colormap gray
        axis equal tight off;
        %caxis([-.5 1]);
        title(['\fontsize{10}' sprintf('diff %d: %.1f kHz', n_fr, stim_type_lookup(n_fr)/1000)]);
    end 
    suptitle('All amps combined');
    savefig(fig1, [ops.data_dir '\Spatial_ave_all ' ops.file_name]);
end


end