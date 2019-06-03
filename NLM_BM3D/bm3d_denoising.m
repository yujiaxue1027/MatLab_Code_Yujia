function returned = bm3d_denoising(inputname,outputname)
t0 = clock;
noisy_data = read_tif_to_mat(inputname);
[rows,cols,layers] = size(noisy_data);
output = zeros(rows,cols,layers);
for i = 1:layers
    [~,tmp] = BM3D(1,noisy_data(:,:,i));
    output(:,:,i) = tmp;
    disp([num2str(i),'/',num2str(layers)]);
end
write_mat_to_tif(uint8(255*output),outputname);
t1 = clock;
delta_t = t1 - t0;
elapsed_time = delta_t(6) + 60*delta_t(5) + 60*60*delta_t(4) + ...
    24*60*60*delta_t(3) + 30*24*60*60*delta_t(2) +...
    365*30*24*60*60*delta_t(1);
disp(['BM3D: Total elapsed time is ',num2str(elapsed_time),' seconds.']);
returned = uint8(255*output);
end

