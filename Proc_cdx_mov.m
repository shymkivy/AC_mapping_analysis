ops.data_dir = 'D:\data\AC\mapping\12_4_21a';
ops.fname = 'AC_mapping1_12_4_21a.cxd'; 

data = bfopen([ops.data_dir '\' ops.fname]);


bin = 2;

[d1, d2] = size(data{1}{1,1});
T = size(data{1},1);

ave_im = zeros(d1,d2);
for n_t = 1:T
    ave_im = ave_im + double(data{1}{n_t,1})/T;
    fprintf('%d-\n', n_t);
end


for n_t = 627:T
    frame = zeros(d1/2, d2/2, 'uint16');
    frame_pre = double(data{1}{n_t,1});
    for n_d1 = 1:d1/2
        d1_source = (n_d1-1)*bin+1;
        for n_d2 = 1:d2/2
            d2_source = (n_d2-1)*bin+1;
            frame(n_d1, n_d2) = mean(mean(frame_pre(d1_source:(d1_source+bin-1), d2_source:(d2_source+bin-1),:),1),2);
        end
    end
    data{1}{n_t,1} = frame;
    fprintf('frame %d\n', n_t);
end

mov2 = zeros(d1/2, d2/2, T, 'uint16');
for n_d1 = 1:d1/2
    for n_d2 = 1:d2/2
        d1_source = (n_d1-1)*bin+1;
        d2_source = (n_d2-1)*bin+1;
        mov2(n_d1, n_d2,:) = mean(mean(mov1(d1_source:(d1_source+bin-1), d2_source:(d2_source+bin-1),:),1),2);
    end
    fprintf('%d/%d', n_d1, d1/2);
end

data{1} = [];

x = data{3};
data{3} = [];

y = data{2};

figure; imagesc(ave_im)

figure; plot(squeeze(mean(mean(mov1,1),2)))

figure; imagesc(data{1}{5000,1})



x2 = x1.dumpXML;

x1.getPlaneCount(0);

x2 = x1(1);

x2(1)





