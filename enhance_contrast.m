function output = enhance_contrast(filename, low_threshold, high_threshold)
t0 = clock;
linear_normalize = @(x) (x - min(x(:)))./(max(x(:))-min(x(:)));
data = read_tif_to_mat(filename);

if isa(data,'uint8')
    data = double(data)./255;
    dtype = 255;
elseif isa(data,'uint16')
    data = double(data)./65535;
    dtype = 65535;
end


if nargin == 1
    low_threshold = 0.01;
    high_threshold = 0.01;
end

if nargin == 2
    high_threshold = low_threshold;
end

hist_obj = histogram(data(:),1024);
bin_edges = hist_obj.BinEdges;
bin_counts = hist_obj.BinCounts;
close all;
pdf = bin_counts ./ sum(bin_counts(:));
cdf = cumsum(pdf);
[~,idx] = min(abs(cdf - low_threshold));
low_clip_val = 0.5 * (bin_edges(idx) + bin_edges(idx + 1));
data(data <= low_clip_val) = low_clip_val;
data = linear_normalize(data);

hist_obj = histogram(data(:),1024);
bin_edges = hist_obj.BinEdges;
bin_counts = hist_obj.BinCounts;
close all;
pdf = bin_counts ./ sum(bin_counts(:));
cdf = cumsum(pdf);
[~,idx] = min(abs(cdf - (1 - high_threshold)));
high_clip_val = 0.5 * (bin_edges(idx) + bin_edges(idx + 1));
data(data >= high_clip_val) = high_clip_val;
output = linear_normalize(data);

output_filename = ['enhanced_',filename];

if dtype == 255
    output = uint8(255*output);
elseif dtype == 65535
    output = uint16(65535*output);
end

write_mat_to_tif(output, output_filename);

t1 = clock;
delta_t = t1 - t0;
elapsed_time = delta_t(6) + 60*delta_t(5) + 60*60*delta_t(4) + ...
    24*60*60*delta_t(3) + 30*24*60*60*delta_t(2) +...
    365*30*24*60*60*delta_t(1);
disp(['Enhance contrast: total elapsed time: ',num2str(elapsed_time),' seconds']);

end

