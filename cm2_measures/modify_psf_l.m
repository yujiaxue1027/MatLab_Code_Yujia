function modified_psf = modify_psf_l(mla_psf,relative_intensity)
% for 2x downsampled psf only
row_range = [1, 300*2, 700*2,972*2];
col_range = [1, 500*2, 900*2, 1296*2];

modified_psf = mla_psf;
for i = 1:3
for j = 1:3
patch = modified_psf(row_range(i):row_range(i+1), col_range(j):col_range(j+1));
patch = patch./sum(patch(:));
modified_psf(row_range(i):row_range(i+1), col_range(j):col_range(j+1)) = patch.*relative_intensity(i,j);
end
end


modified_psf = modified_psf./sum(modified_psf(:));

end
