function [est_bg,est_fg] = rpca_decomposition(filename,lambda,combined, start_frame, end_frame)
t0 = clock;
data = read_tif_to_mat(filename);

if isa(data,'uint8')
    data = double(data)./255;
    dtype = 255;
elseif isa(data,'uint16')
    data = double(data)./65535;
    dtype = 65535;
end

[~,~,f] = size(data);

if nargin < 2
    lambda = 0.001;
end

if nargin < 3
    combined = 2;
end

if nargin < 4
    start_frame = 1;
    end_frame = f;
end

data = data(:,:,start_frame:end_frame);
[r,c,f] = size(data);

data_flatten = reshape(data,[r*c, f]);
[est_lr, est_sp, ~] = inexact_alm_rpca(data_flatten, lambda);
est_bg = reshape(est_lr,[r,c,f]);
est_fg = reshape(est_sp,[r,c,f]);

bg_filename = ['est_bg_',num2str(lambda),'_',filename];
fg_filename = ['est_fg_',num2str(lambda),'_',filename];
combined_filename = ['est_',num2str(lambda),'_',filename];

if dtype == 255
    est_bg = uint8(255*est_bg);
    est_fg = uint8(255*est_fg);
elseif dtype == 65535
    est_bg = uint16(65535*est_bg);
    est_fg = uint16(65535*est_fg);
end

if combined == 0
    write_mat_to_tif(est_bg, bg_filename);
    write_mat_to_tif(est_fg, fg_filename);
elseif combined == 1
    write_mat_to_tif(cat(1,est_bg, est_fg), combined_filename);
elseif combined == 2
    write_mat_to_tif(cat(2,est_bg, est_fg), combined_filename);
end

t1 = clock;
delta_t = t1 - t0;
elapsed_time = delta_t(6) + 60*delta_t(5) + 60*60*delta_t(4) + ...
    24*60*60*delta_t(3) + 30*24*60*60*delta_t(2) +...
    365*30*24*60*60*delta_t(1);
disp(['Robust-PCA: total elapsed time: ',num2str(elapsed_time),' seconds']);

end

