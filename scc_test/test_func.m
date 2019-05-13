function output = test_func(rows,cols, save_path_and_name)
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));

x = rand([rows,cols]);
fx = F(x);
output = fx;
save(save_path_and_name, 'fx');

end

