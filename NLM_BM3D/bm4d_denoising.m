function returned = bm4d_denoising(inputname,outputname)
t0 = clock;
noisy_data = read_tif_to_mat(inputname);
output = bm4d(noisy_data, 'Gauss');
write_mat_to_tif(uint8(output),outputname);
t1 = clock;
delta_t = t1 - t0;
elapsed_time = delta_t(6) + 60*delta_t(5) + 60*60*delta_t(4) + ...
    24*60*60*delta_t(3) + 30*24*60*60*delta_t(2) +...
    365*30*24*60*60*delta_t(1);
disp(['BM4D: Total elapsed time is ',num2str(elapsed_time),' seconds.']);
returned = uint8(output);
end

