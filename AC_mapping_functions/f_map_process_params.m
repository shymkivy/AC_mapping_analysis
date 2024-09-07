function params = f_map_process_params(ops)

%kernel = ones(kernel_size)/kernel_size^2;
kernel = hamming(ops.kernel_size)*hamming(ops.kernel_size)';
params.kernel = kernel/sum(kernel(:));

% figure;
% surf(kernel)

%% Load stim parameters
stim_params = load([ops.data_dir '\' ops.file_name,'_stim_data.mat']);
params.stim_params = stim_params;

params.freq_amp_lookup = reshape((1:stim_params.num_freqs*numel(stim_params.modulation_amp))',stim_params.num_freqs, numel(stim_params.modulation_amp));

%% compute used frequencies
params.num_trial_types = stim_params.num_freqs*numel(stim_params.modulation_amp);

if isfield(stim_params, 'all_trial_types')
    stim_type_lookup = stim_params.all_trial_types;
else
    control_carrier_freq = zeros(stim_params.num_freqs, 1);
    control_carrier_freq(1) = stim_params.start_freq;
    for n_fr = 2:stim_params.num_freqs
        control_carrier_freq(n_fr) = control_carrier_freq(n_fr-1) * stim_params.increase_factor;
    end

    stim_type_lookup = zeros(params.num_trial_types,numel(stim_params.modulation_amp),2);
    for n_amp = 1:numel(stim_params.modulation_amp)
        stim_type_lookup(:,n_amp,1) = control_carrier_freq;
        stim_type_lookup(:,n_amp,2) = stim_params.modulation_amp(n_amp);
    end
    stim_type_lookup = reshape(stim_type_lookup, [], 2);
end
params.stim_type_lookup = stim_type_lookup;

%% compute what window to use based on isi and stim duration

baseline_frames = round(ceil((-ops.baseline_window_time(1))*ops.framerate));
sta_num_frames = round(baseline_frames+ceil((stim_params.duration + stim_params.isi)*ops.framerate));

params.alignment_mode = ops.alignment_mode;
if ops.alignment_mode == 1
    params.prepulse_frame_shift = (stim_params.pretrial_pulse_duration + stim_params.pretrial_pulse_wait) * ops.framerate;
end


if stim_params.pretrial_LED_pulse
    if ops.show_alignment_pulse
        sta_num_frames = sta_num_frames + (stim_params.pretrial_pulse_duration + stim_params.pretrial_pulse_wait) * framerate;
    else
        % remove some time to avoid the pretrial pulse in analysis 
        sta_num_frames = sta_num_frames - (0.3) * framerate;
    end
end


% signal starts in frame 0 of sta_window
sta_window = (0:sta_num_frames-1)-baseline_frames;
% frame 0 of sta_window corresponds to first time point after t=0
t_sta = (sta_window+1)/ops.framerate; % in sec

baseline_window_frames = logical((t_sta>ops.baseline_window_time(1)) .* (t_sta<ops.baseline_window_time(2)));

onset_window_frames = logical((t_sta>ops.onset_window_time(1)) .* (t_sta<ops.onset_window_time(2)));

% offset_window_frames = logical((time_stim_window>offset_window_time(1)) .* (time_stim_window<offset_window_time(2)));

%%
params.sta_window = sta_window;
params.t_sta = t_sta;
params.baseline_window_frames = baseline_window_frames;
params.onset_window_frames = onset_window_frames;

end