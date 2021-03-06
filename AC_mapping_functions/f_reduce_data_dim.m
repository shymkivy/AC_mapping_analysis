function f_reduce_data_dim(data_out, params)

[d1,d2,t_bin,trials] = size(data_out.vec_frame_data);
data1 = double(reshape(data_out.vec_frame_data(:,:,:,data_out.trial_types == 4),d1*d2,[]));
data_mean = mean(data1,2);
data1_norm = (data1 - data_mean)./std(data1,[],2);
[U,S,V] = svd(data1_norm);

[W,H] = nnmf(data1_norm,15);

%% check which comps are modified by signal
base_frames = params.t_sta<=0;

for n_comp = 1:15
    trial_V = reshape(H(n_comp,:),t_bin,[]);
    base_sig = trial_V(base_frames,:);
    trial_base_mean = mean(base_sig);
    base_sig_norm = (base_sig - trial_base_mean);
    trial_base_std = std(base_sig_norm(:));
    figure; 
    subplot(2,1,1)
    imagesc(reshape(W(:,n_comp),d1,d2));
    title(['comp ' num2str(n_comp)]);
    subplot(2,1,2)
    hold on; axis tight;
    plot(params.t_sta, trial_V-trial_base_mean, 'color', [.6 .6 .6])
    plot(params.t_sta, mean(trial_V-trial_base_mean,2), 'm', 'LineWidth',2)
    plot(params.t_sta, ones(1,t_bin)*3*trial_base_std, '--r')
    plot(params.t_sta, ones(1,t_bin)*-3*trial_base_std, '--r')
   
end


      

%% plot
n_comps = 6;
data_LR = U(:,n_comps)*S(n_comps,n_comps)*V(:,n_comps)';

vid_stack = f_make_single_trial_tiles(trial_vec_data, params.onset_window_frames);
        

vid_stack = reshape(data_LR,d1,d2,[]);
f_save_tif_stack2_YS(vid_stack, 'test.tif');

figure; imagesc(reshape(U(:,n_comp),d1,d2))

n_comp = 14;
trial_V = reshape(V(:,n_comp),t_bin,[]);
trial_base = mean(trial_V(1:9,:));
figure; hold on; axis tight;
plot(params.t_sta, trial_V-trial_base, 'color', [.6 .6 .6])
plot(params.t_sta, mean(trial_V-trial_base,2), 'm', 'LineWidth',2)
title(['comp ' num2str(n_comp)])

figure; plot(diag(S))


end