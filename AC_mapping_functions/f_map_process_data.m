function data_out = f_map_process_data(ops, params)

% load data or preprocess
if exist([ops.data_dir '\' ops.file_name,'_vectorized_data_align2.mat'], 'file') == 2
    disp('Loading preprocessed dataset _vectorized_data_align2...');
    data_out = load([ops.data_dir '\' ops.file_name,'_vectorized_data_align2.mat']);
elseif exist([ops.data_dir '\' ops.file_name,'_vectorized_data_align1.mat'], 'file') == 2
    disp('Loading preprocessed dataset _vectorized_data_align1...');
    data_out = load([ops.data_dir '\' ops.file_name,'_vectorized_data_align1.mat']);
else

    %% First load movie and get ave trace of all pixels
    % load alignment data
    if exist([ops.data_dir '\' ops.file_name,'_processing_params.mat'], 'file')
        proc_params = load([ops.data_dir '\' ops.file_name,'_processing_params.mat']);
    else
        save([ops.data_dir '\' ops.file_name,'_processing_params.mat'], '');
        proc_params = struct;
    end
    
    
    disp('Loading video...');
    data_Y=bigread([ops.data_dir '\' ops.file_name '.tif']);
    disp('Finished loading the video');
    
    [d1, d2, dT] = size(data_Y);
    
    if isfield(proc_params, 'FOV')
        FOV = proc_params.FOV;
    else
        FOV = struct;
    end
    
    if ~isfield(FOV, 'mean_Y_raw')
        FOV.mean_Y_raw = mean(data_Y,3);
    end
    
    if ops.reduce_pix_num
        if ~isfield(FOV, 'xPix')
            figure;
            imagesc(FOV.mean_Y_raw);
            axis equal tight;
            title('Mean sig full, Select active region (2 clicks)');
            [FOV.xPix, FOV.yPix] = ginput(2);
            title('Mean sig full');         
            FOV.xPix(1) = max(FOV.xPix(1),1);
            FOV.xPix(2) = min(FOV.xPix(2),d2);
            FOV.yPix(1) = max(FOV.yPix(1),1);
            FOV.yPix(2) = min(FOV.yPix(2),d1);
            FOV.xPix = round(FOV.xPix);
            FOV.yPix = round(FOV.yPix);
        end
        data_Y = data_Y(FOV.yPix(1):FOV.yPix(2),FOV.xPix(1):FOV.xPix(2),:);
    end
    FOV.mean_Y_cut = mean(data_Y,3);
    figure;
    imagesc(FOV.mean_Y_cut);
    axis equal tight;
    title('Mean sig cut');
    
    proc_params.FOV = FOV;
    if ~exist([ops.data_dir '\' ops.file_name,'_processing_params.mat'], 'file')
        save([ops.data_dir '\' ops.file_name,'_processing_params.mat'],'FOV');
    else
        save([ops.data_dir '\' ops.file_name,'_processing_params.mat'],'FOV', '-append');
    end
    
    [d1, d2, dT] = size(data_Y);
    %% Voltage 
    % Load and process voltage traces
    data_voltage = csvread([ops.data_dir '\' ops.file_name '.csv'], 1, 0);

    % process DAQ voltage for stim type channel
    volt_data_filt = zeros(size(data_voltage,1),3);
    volt_data_filt(:,1) = medfilt1(data_voltage(:,ops.freq_volt_ch+1), 49);
    volt_data_filt(:,2) = medfilt1(data_voltage(:,ops.LED_ch+1), 49)/3;
    volt_data_filt(:,3) = data_voltage(:,ops.locomotion_ch+1);
    figure;
    labels1 = {'Stim ch', 'LED ch', 'Loco ch'};
    for n_pl = 1:3
        subplot(3,1,n_pl);
        plot(data_voltage(:,1)/1000,volt_data_filt(:,n_pl));
        title(labels1{n_pl})
    end

%% Alignments
    proc_params.alignment_mode = ops.alignment_mode;
    
    if ops.alignment_mode == 2
        disp('Alignment mode 2, using camera timing data...');
        %% Use frame time data from camera to fill missing frames

        % load CCD data
        data_CCD = csvread([ops.data_dir '\' ops.file_name '_CCD_out.csv'], 1, 0);
        t_frames_CCD = data_CCD(:,3);
        clear data_CCD;

        % now fill in missing frames
        disp('Interpolating video to fill missing frames...');

        % find skipped frames first
        frame_times = [0; diff(t_frames_CCD)];
        %frame_period_appx = mode(frame_times);      % it should be the most common frame time
        frame_period_appx = mode(medfilt1(diff(t_frames_CCD),15));
        skipped_frames = [];
        for ii = 1:numel(frame_times)
            if frame_times(ii) > 1.5*frame_period_appx
                skipped_frames = [skipped_frames; ii, round(frame_times(ii) / frame_period_appx - 1)];
            end
        end
        if sum(skipped_frames(:)) > 0
            figure; plot(diff(t_frames_CCD));
            title('Video frame times');
            [data_Y_full, t_frames] = if_interp_skipped_frames(data_Y, skipped_frames, t_frames_CCD);
        else
            disp('No skipped frames')
            t_frames = t_frames_CCD;
            data_Y_full = single(data_Y);
        end
        clear data_Y;
        
        params.framerate = 1000*mode(diff(t_frames));
        
        % check if fill was good 
        figure;
        plot(diff(t_frames));
        title('Check frame period after filling frames');
        ylabel('frame period');
        xlabel('frame number')

        clear t_frames_CCD;

        %% normalize mean ca trace
        % normalize mean Ca signal from frames
        mean_frames = squeeze(mean(mean(data_Y_full,1),2));
        base_sub_mean_frames = mean_frames - min(mean_frames);
        norm_mean_Ca_data = base_sub_mean_frames/(max(base_sub_mean_frames));
        clear mean_frames base_sub_mean_frames;
        %save([file_name '_mean_Ca_trace'],'norm_Ca_data_filled');

        % extract other parameters
        num_frames = numel(t_frames);
        frame_period = mean(diff(t_frames));



        %% extract times of alignment pulses

        disp('Alignment...')

        if ~isfield(proc_params, 'synch_pulse_frames')

            figure;
            plot(norm_mean_Ca_data); axis tight;
            title('Click before and after alignment pulses (4 clicks)');
            legend('Mean Ca trace from video');
            [x,~]=ginput(4);
            close;
            
            synch_pulse_frames = if_extract_artifact(norm_mean_Ca_data, round(x));

            proc_params.synch_pulse_frames = synch_pulse_frames;
            save([ops.data_dir '\' ops.file_name,'_processing_params.mat'],'synch_pulse_frames', '-append');
        end

        % exract artifact locations in LED voltage trace
        if ~isfield(proc_params, 'synch_pulse_ms')

            figure;
            plot(volt_data_filt(:,2)); axis tight;
            title('Click before and after alignment artifacts (4 clicks)');
            legend('Voltage, LED channel');
            [x,~]=ginput(4);
            close;

            % extract trigger location in voltage trace
            synch_pulse_ms = if_extract_artifact(volt_data_filt(:,2), round(x));
            proc_params.synch_pulse_ms = synch_pulse_ms;
            save([ops.data_dir '\' ops.file_name,'_processing_params.mat'],'synch_pulse_ms', '-append');
        end


    %     % get alignment paramteters
    %     try
    %         load(strcat(file_name,'_alignment2'),'synch_coord');
    %     catch
    %         % plot traces before fixing the shift
    %         figure;
    %         dcm_obj = datacursormode;
    %         set(dcm_obj,'UpdateFcn', @NewCallback_YS)
    %         plot(t_ms,daq_freq_voltage/max(daq_freq_voltage));
    %         hold on;
    %         plot(t_frames,norm_Ca_data_filled);
    %         title('Traces before correction');
    %         legend('DAQ voltage trace', 'Mean Ca trace from video');
    % 
    %         synch_coord = zeros(2,2);
    % 
    %         synch_coord(1,1) = input('Input first synch peak start of Ca trace (last point before rise):');
    %         synch_coord(2,1) = input('Input first synch peak start of voltage trace (last point before rise):');
    %         synch_coord(1,2) = input('Input second synch peak start of Ca trace (last point before rise):');
    %         synch_coord(2,2) = input('Input second synch peak start of voltage trace (last point before rise):');
    %         close;
    %         % frames_voltage_peak_start = [5400 492200; 10128 511049];
    % 
    %         save(strcat(file_name,'_alignment2'),'synch_coord');
    %     end


        %% align

        % if exist('alignment', 'var')
        %     fprintf('Loaded saved alignment data\n');
        % else
        synch_coord = [t_frames(proc_params.synch_pulse_frames)*1000,proc_params.synch_pulse_ms]';

        alignment.scaling_factor = (synch_coord(1,2)-synch_coord(1,1))/(synch_coord(2,2)-synch_coord(2,1));
        alignment.shift = round(synch_coord(1,1) - (synch_coord(2,1)*alignment.scaling_factor));
        
        
        % save([file_out '_alignment2'], 'alignment', '-append');
        % fprintf('Saved alignment data\n');

        % align here
        volt_data_aligned = align_volt_by_scale_shift(volt_data_filt, alignment.scaling_factor, alignment.shift);

        % trim or fill voltage to make make length the same with Ca video
        if round(1000*(t_frames(end) + frame_period) + 1) > size(volt_data_aligned,1)
            temp_fill = zeros(round(1000*(t_frames(end) + frame_period) + 1 - size(volt_data_aligned,1)), size(volt_data_aligned,2));
            volt_data_aligned = [volt_data_aligned; temp_fill];
            clear temp_fill;
        else
            volt_data_aligned = volt_data_aligned(1:round(1000*(t_frames(end) + frame_period) + 1),:);
        end
        clear volt_data_filt;


        %% bin voltage data into frame bins

        volt_data_binned = zeros(num_frames, 3);

        for frame=1:num_frames
            frame_start_index = round(1000*t_frames(frame))+1;
            frame_end_index = round(1000*(t_frames(frame)+frame_period));

            % averages visual stim voltage trigger
            volt_data_binned(frame, 1)=median(volt_data_aligned(frame_start_index:frame_end_index,1),1);
            
            % averages LED voltage trigger 
            volt_data_binned(frame, 2)=median(volt_data_aligned(frame_start_index:frame_end_index,2),1);
            % takes the absolute val of first derivative
            volt_data_binned(frame, 3)=mean(volt_data_aligned(frame_start_index:frame_end_index,3),1);

        end
        volt_data_binned(:, 1) = round(volt_data_binned(:, 1)/4*params.num_trial_types);
        % process the movement channel
        volt_data_binned(:, 3) = abs(gradient(volt_data_binned(:, 3)));
        clear frame frame_start_index frame_end_index volt_data_aligned;


        %% Get mapping phase
        disp('Extracting stim timing...');
        % select mapping phase 
        if ~isfield(proc_params, 'mapping_phase_x')
            figure;
            hold on;
            plot(volt_data_binned(:, 1));
            plot(norm_mean_Ca_data); axis tight;
            title('Pick start and end of mapping of phase.');
            [mapping_phase_x,~]=ginput(2);
            close;
            proc_params.mapping_phase_x = mapping_phase_x;
            save([ops.data_dir '\' ops.file_name,'_processing_params.mat'],'mapping_phase_x', '-append');
        end
        
        mapping_phase = zeros(size(volt_data_binned,1),1);
        mapping_phase(round(proc_params.mapping_phase_x(1)):round(proc_params.mapping_phase_x(2))) = 1;

        %% Extract stim times

        stim_frame = zeros(num_frames,1);
        stim_cutoff = 1/2;
        for frame = 2:num_frames
            if mapping_phase(frame) == 1
                if (volt_data_binned(frame) > stim_cutoff) && (volt_data_binned(frame-1) < stim_cutoff)
                    stim_frame(frame) = 1;
                end
            end
        end
        clear stim_cutoff;
        % extract stim times
        stim_frame_index = find(stim_frame == 1);
        num_trials = sum(stim_frame);
        fprintf('%d trials extracted\n', num_trials);
        clear stim_frame;

        %% Extract sta data
        disp('Extracting data...');
        % Extract stimulus direction
        trial_types = zeros(num_trials,1);
        for trial = 1:num_trials
            % do a median over frames to prevent some unlikely noise
            trial_types(trial) = median(volt_data_binned(stim_frame_index(trial):stim_frame_index(trial)+round(params.stim_params.duration/frame_period)));
        end
        %trial_type_orientations = rem(trial_type_directions-1, 6)+1;


        % specifies which stim type was given trial
        vec_frame_data = uint16(zeros([d1, d2, numel(params.sta_window), params.num_trial_types]));

    %     loco_trials = zeros(num_trials,1);

        for trial = 1:num_trials
            % copies the Ca signal of given window
            temp_frame = stim_frame_index(trial);
            vec_frame_data(:,:,:,trial) = data_Y_full(:,:, temp_frame+params.sta_window);

            % mark locomotion trials
    %         if sum(volt_data((temp_frame-baseline_window_frames):(temp_frame+stim_ave_window_frames-1),3)) > 0
    %             loco_trials(trial) = 1;
    %         end
        end
        %%

        %check alignment
        figure;
        plot(t_frames, norm_mean_Ca_data);
        hold on;
        plot(t_frames, volt_data_binned(:,1));
        plot(t_frames, volt_data_binned(:,2));
        plot(t_frames, mapping_phase*10);  axis tight;
        title('Check if alignment is good');
        legend('Ca trace', 'Binned, aligned, indexed voltage', 'LED', 'Mapping phase');
        
        data_out.vec_frame_data = vec_frame_data;
        data_out.trial_types = trial_types;
        data_out.volt_data_binned = volt_data_binned;
        disp('Saving data output...');
        save([ops.data_dir, '\' ops.file_name '_vectorized_data_align2'], 'vec_frame_data', 'trial_types', 'volt_data_binned', '-v7.3');
    end

%% Alignment 1

    if ops.alignment_mode == 1
        % first specify the mapping period
        disp('Alignment mode 1, using light artifacts before each stim...');

        t_frames = (1:1:dT)*1000/framerate;
        %%
        disp('Extracting artifact triggered average...');

        if exist([data_dir '\' file_name,'_alignment1.mat'], 'file') == 2
            load([data_dir '\' file_name,'_alignment1.mat']);
        end

        if ~exist('mapping_phase_bounds_frames', 'var')
            figure;
            plot(norm_mean_Ca_data);
            title('Pick start and end of mapping of phase');
            legend('Mean Ca trace from video');
            [mapping_phase_bounds_frames,~]=ginput(2);
            close;
            mapping_phase_bounds_frames = round(mapping_phase_bounds_frames);
            save([data_dir '\' file_name,'_alignment1.mat'],'mapping_phase_bounds_frames');
        end

        % provide experimetnal phase
        mapping_phase_frames = zeros(size(norm_mean_Ca_data));
        mapping_phase_frames(mapping_phase_bounds_frames(1):mapping_phase_bounds_frames(2)) = 1;

        figure;
        plot(norm_mean_Ca_data)
        hold on;
        plot(mapping_phase_frames);



        %% artifact triggered average


        % define the mapping phase
        if ~exist('mapping_phase_bounds_volt', 'var')
            figure;
            plot(volt_data_filt(:,2));
            title('Pick start and end of mapping of phase.');
            [mapping_phase_bounds_volt,~]=ginput(2);
            close;
            volt_mapping_stimcuts = round(mapping_phase_bounds_volt);
            save([data_dir '\' file_name,'_alignment1.mat'],'mapping_phase_bounds_volt', '-append');
        end

        % here we are extracting the trial types
        mapping_phase_volt = zeros(size(volt_data_filt(:,2)));
        mapping_phase_volt(volt_mapping_stimcuts(1):volt_mapping_stimcuts(2)) = 1;


        %% here extract the trial types    
        daq_freq_voltage_stim_times = zeros(size(indexed_daq_freq_voltage));
        for ii = 2:length(indexed_daq_freq_voltage)
            if mapping_phase_volt(ii) > 0
                if (indexed_daq_freq_voltage(ii) > 0.5) && (indexed_daq_freq_voltage(ii-1) < 0.5)
                    daq_freq_voltage_stim_times(ii) = max(indexed_daq_freq_voltage(ii:round(ii+1000*duration)));
                end
            end
        end
        fprintf('Total trial extracted from voltage %d \n', sum(sign(daq_freq_voltage_stim_times)));
        trial_types = daq_freq_voltage_stim_times(daq_freq_voltage_stim_times > 0);



        % pull out atrifact times
        atrif_index_frames = zeros(size(norm_mean_Ca_data));

        if ops.inverted_alignment_pulses == 1
            % flip the signal to fit the original code
            artifact_processing_sig = -norm_mean_Ca_data;
        else
            artifact_processing_sig = norm_mean_Ca_data;
        end

        for ii = 2:size(norm_mean_Ca_data,1)
            if mapping_phase_frames(ii) > 0
                if (artifact_processing_sig(ii) - artifact_processing_sig(ii-1) > 0.1) && (artifact_processing_sig(ii-1) - artifact_processing_sig(ii-2) < 0.1)
                    atrif_index_frames(ii) = 1;
                end
            end
        end
        % pull out artifact times index
        artif_index_times = find(atrif_index_frames == 1);
        % get the frame when stim started
        stim_frame_index = artif_index_times + prepulse_frame_shift;
        num_trials = sum(atrif_index_frames);
        fprintf('Total artifacts extracted %d \n', num_trials);


        figure;
        plot(norm_mean_Ca_data);
        hold on;
        plot(mapping_phase_frames);
        plot(atrif_index_frames/2);
        title(sprintf('Atrifact triggers, %d total', num_trials));
        legend('CCD averaged signal', 'Experimental mapping phase', 'Pulled out artifact');


        % extracting data from movie
        % specifies which stim type was given trial
        vec_frame_data = uint16(zeros([d1, d2, sta_num_frames, num_trials]));

        for trial = 1:length(stim_frame_index)
            % copies the Ca signal of given window
            temp_frame = stim_frame_index(trial);
            vec_frame_data(:,:,:,trial) = data_Y_full(:,:, temp_frame+sta_window);
        end
        

        % throw out trials with lost frames, when camera fucks up
        for remove_frames = 1
            if remove_frames == 1
                lower_removal_threshold = 1;
                upper_removal_threshold = 3;

                figure;
                plot(diff(artif_index_times));
                hold on;
                plot(ones(1,length(diff(artif_index_times)))*(median(diff(artif_index_times))+upper_removal_threshold), 'r');
                plot(ones(1,length(diff(artif_index_times)))*(median(diff(artif_index_times))-lower_removal_threshold), 'r');
                legend('Frame number', 'exclusion thresholds');
                title('Frame number between light artifacts');
                xlabel('Trial');
                ylabel('Frame number');

                frames_to_remove = [or((diff(artif_index_times) <= (median(diff(artif_index_times))-lower_removal_threshold)),((diff(artif_index_times) >= (median(diff(artif_index_times))+upper_removal_threshold)))); 1];
                fprintf('%d frames removed from analysis\n', sum(frames_to_remove));

                % remove frames
                trial_types(logical(frames_to_remove)) = [];
                vec_frame_data(:,:,:,logical(frames_to_remove)) = [];
            end
        end


        % save output
        toc
        disp('Saving data output...');
        save([data_dir '\' file_name '_vectorized_data_align1'],'vec_frame_data', 'vec_ave_data', 'trial_types', '-v7.3');
        %clear data_raw;
    end
end


end

%%
function [data_Y_full, t_frames] = if_interp_skipped_frames(data_Y, skipped_frames, t_frames_CCD);

num_skipped_frames = sum(skipped_frames(:,2));
fprintf('%d frames missing\n', num_skipped_frames);

% fill in skipped frames
t_frames = t_frames_CCD;
temp_shift = 0;
for ii = 1:size(skipped_frames,1)
    current_frame = skipped_frames(ii,1) + temp_shift;

    temp_fill_time = linspace(t_frames(current_frame-1), t_frames(current_frame), skipped_frames(ii,2) + 2)';
    t_frames = [t_frames(1:current_frame-1);temp_fill_time(2:end-1); t_frames(current_frame:end)];

    temp_shift = temp_shift+ skipped_frames(ii,2);
end

% interpolate Ca videos along T axis to fill frames
disp('Interpolating...')
X = 1:size(data_Y,2);
Y = (1:size(data_Y,1))';
T = t_frames_CCD;
Tq = t_frames;

data_Y_full = interp3(X,Y,T,single(data_Y),X,Y,Tq, 'linear');
clear X Y T Tq V;

figure;
hold on;
plot(t_frames, squeeze(mean(mean(data_Y_full))), '-o');
plot(t_frames_CCD, squeeze(mean(mean(data_Y))), 'o-');
title(sprintf('Check the interpolation; %d skipped frames\n', num_skipped_frames));
legend('After', 'Before');
xlabel('Time, sec');
ylabel('Mean signal');

end

function synch_pulse_frames = if_extract_artifact(trace_in, coord)
    
    synch_pulse_positive = mean(trace_in) < 0.5;

    num_pulses = round(numel(coord)/2);
    synch_pulse_frames = zeros(num_pulses,1);
    
    % extract pulse rise coordinates in trace
    for ii = 1:num_pulses
        pulse_time = zeros(size(trace_in));
        pulse_time(coord(ii*2-1):coord(ii*2)) = 1;
        
        diff_trace_in = [0; diff(trace_in)].*pulse_time;
        
        if synch_pulse_positive
            synch_pulse_frames(ii) = find(diff_trace_in == max(diff_trace_in));
        else
            synch_pulse_frames(ii) = find(diff_trace_in == min(diff_trace_in));
        end
    end
    
end
