function f_map_plot_trial_ave(vec_frame_data, trial_types, params, ops)
%% get mean ca data across all pixels
vec_ave_data = squeeze(mean(mean(vec_frame_data,1),2));
vec_ave_data = vec_ave_data - min(vec_ave_data(:));
vec_ave_data = vec_ave_data / max(vec_ave_data(:));

%% select trials
num_frames = size(vec_frame_data,3);
freq_amp_lookup = params.freq_amp_lookup;
modulation_amp = params.stim_params.modulation_amp;
num_freqs = params.stim_params.num_freqs;
duration = params.stim_params.duration;
stim_type_lookup = params.stim_type_lookup;

%% compute averages over all trials
%analysis_window = analysis_window - 5
traces_ave = zeros(num_frames, params.num_trial_types);
traces_sem = zeros(num_frames, params.num_trial_types);

sp_dim = [ceil(num_freqs/5) 5];

% define windows
for n_tr = 1:params.num_trial_types
    temp_trials = vec_ave_data(:,trial_types == n_tr);
    % baseline subtraceted average
    traces_ave(:,n_tr) = mean(temp_trials,2) - mean(mean(temp_trials(params.baseline_window_frames, :),1));
    traces_sem(:,n_tr) = std(temp_trials,[],2)/sqrt(sum(trial_types == n_tr)-1);
end
y_min = min([traces_ave(:)-traces_sem(:); 0]);
y_max = max(traces_ave(:)+traces_sem(:));

nfr_colors = {'k', 'b', 'g', 'm'};
sp = cell(numel(modulation_amp),1);
legend_tag = cell(numel(modulation_amp),1);
for n_amp = 1:numel(modulation_amp)
    legend_tag{n_amp} = ['amp ' num2str(modulation_amp(n_amp))];
end
figure;
for n_fr = 1:num_freqs
    subplot(sp_dim(1),sp_dim(2),n_fr); hold on;
    for n_amp = 1:numel(modulation_amp)
        temp_ave_data = traces_ave(:,freq_amp_lookup(n_fr,n_amp));
        temp_sem_data = traces_sem(:,freq_amp_lookup(n_fr,n_amp));
        sp{n_amp} = shadedErrorBar(params.t_sta, temp_ave_data, temp_sem_data,'lineprops',nfr_colors{n_amp});      
    end
    patch([0 0 duration duration],[y_min y_max y_max y_min] ,[1 .6 .6], 'LineStyle', 'none', 'FaceAlpha', 0.3);
    xlabel('Time, sec');
    axis tight;
    title(['\fontsize{10}' sprintf('%d: %.1f kHz', n_fr, stim_type_lookup(n_fr)/1000)]);
    if n_fr == num_freqs
        leg_patch = [sp{:}];
        legend([leg_patch.patch], legend_tag,'Location','southeast')
    end
end
suptitle('Mean sig across freqs and amps');

%% avarage over all amps

if numel(modulation_amp)>1
    traces_ave = zeros(num_frames, num_freqs);
    traces_sem = zeros(num_frames, num_freqs);

    % define windows
    for n_fr = 1:num_freqs 
        temp_trials = vec_ave_data(:,logical(sum(trial_types == freq_amp_lookup(n_fr,:),2)));
        % baseline subtraceted average
        traces_ave(:,n_fr) = mean(temp_trials,2) - mean(mean(temp_trials(params.baseline_window_frames, :),1));
        traces_sem(:,n_fr) = std(temp_trials,[],2)/sqrt(sum(sum(trial_types == freq_amp_lookup(n_fr,:),2))-1);
    end
    y_min = min([traces_ave(:)-traces_sem(:); 0]);
    y_max = max(traces_ave(:)+traces_sem(:));

    figure;
    for n_fr = 1:num_freqs
        subplot(sp_dim(1), sp_dim(2),n_fr); hold on;
        sp{1} = shadedErrorBar(params.t_sta, traces_ave(:,n_fr), traces_sem(:,n_fr),'lineprops','k');      
        patch([0 0 duration duration],[y_min y_max y_max y_min] ,[1 .6 .6], 'LineStyle', 'none', 'FaceAlpha', 0.3);
        xlabel('Time, sec');
        axis tight;
        title(['\fontsize{10}' sprintf('%d: %.1f kHz', n_fr, stim_type_lookup(n_fr)/1000)]);
    end
    suptitle('Mean sig over all amps');
end
end