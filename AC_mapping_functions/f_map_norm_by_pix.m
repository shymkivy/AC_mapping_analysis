function data_out = f_map_norm_by_pix(data)

[d1,d2,~,~] = size(data);

data = single(reshape(data, d1 * d2, []));

% first normalize each pizel by mean signal
data = data ./ mean(data,2);

% now do standard normalization
data = data - min(data(:));

data = data / max(data(:));

data_out = single(reshape(data, d1, d2, []));

end