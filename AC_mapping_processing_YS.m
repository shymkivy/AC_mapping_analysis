%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%   A1 mapping: this processes mapping data 
%
%   Input:  CCD movie tiff, preferably 256x256 or fewer pixels 
%           _stim_data.mat file with info on stimuli
%           
%           Alignment mode 1:
%               pretrial pulses duting stim presentation
%
%           Alignment mode 2:
%               _CCD_out.csv file with frame times from Hamamatsu software
%
%   Output: _vectorized_data_align1/2 - sorted data by trials
%
%   Processing steps:
%       blah blah
%
%   Required functions(might be internal version at end): 
%       synch_pulse_frames
%       align_volt_trace
%
%   Last update: Yuriy 5/6/20
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
close all;

%% Processing parameters
addpath('C:\Users\rylab_dataPC\Desktop\Yuriy\AC_mapping')
addpath('C:\Users\rylab_dataPC\Desktop\Yuriy\AC_mapping\AC_mapping_functions')

ops.data_dir = 'F:\data\Auditory\2018\10_2_18_mapping';
%ops.data_dir = 'L:\data\Auditory\2018\11_20_18_mapping';

ops.file_name = 'AC_mapping5_9_28_18';

%file_name = 'A1_freq_sweeps2_9_28_18';

%parameters_file_name = 'A1_mapping2_9_20_18';
ops.save_trial_analysis = 0; % not ready
ops.save_single_trial_vids = 1;
ops.save_trial_ave_vids = 1;

% select alignment mode
% = 1 for artifact triggered average;
% = 2 to use with CCD csv output
ops.alignment_mode = 2; % 1 or 2
ops.reduce_pix_num = 1;

% what type of excitation, xenon lamp with LED artifact, or LED excitation
% alone? LED alone turns off to generate light artifact
ops.inverted_alignment_pulses = 1; % = 0; for xenon lamp with LED

% Specify channels
ops.freq_volt_ch = 1;
ops.LED_ch = 2;
ops.locomotion_ch = 3;



%% Analysis parameters
ops.show_alignment_pulse = 0; % in the final plots, do you want to show it
ops.framerate = 30;                     % frames per second

ops.baseline_window_time = [-.3 -0.005];        % in sec
ops.onset_window_time = [.05 .450];             % in sec
%offset_window_time = [550 950]/1000;        % in sec

% Analysis parameters
ops.smooth_frames = 0; % % during analysis, smoothing frames via convolution
ops.kernel_size = 30; % smoothing parameters

ops.subplot_dimensions = [2 5];

%%
params = f_map_process_params(ops);

%% Align and extract data
data_out = f_map_process_data(ops, params);

%% Analysis
f_map_analysis(data_out, params, ops);

%%
%f_reduce_data_dim(data_out, params);
%%

disp('Done');
