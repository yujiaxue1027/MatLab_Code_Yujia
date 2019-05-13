function [row_co, col_co] = locate_mass_center(patch)
VEC = @(x) x(:);
[rows,cols] = size(patch);
[col_grid, row_grid] = meshgrid([1:1:cols],[1:1:rows]);
row_co = round(sum(VEC(row_grid.*patch))./sum(VEC(patch)));
col_co = round(sum(VEC(col_grid.*patch))./sum(VEC(patch)));
end

