function f_reduce_data_dim(data_out, params)

[d1,d2,t_bin,trials] = size(data_out.vec_frame_data);
data1 = double(reshape(data_out.vec_frame_data(:,:,:,data_out.trial_types == 2),d1*d2,[]));
data_mean = mean(data1,2);
data1_norm = (data1 - data_mean)./std(data1,[],2);
[U,S,V] = svd(data1_norm);

n_comps = 6;
data_LR = U(:,n_comps)*S(n_comps,n_comps)*V(:,n_comps)';

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