function output = M_dagger(input,M,threshold)
numerator = sum(input.*M, 3);
denominator = sum(M.^2, 3) + eps;
tmp_mask = denominator <= threshold; %% mask = 1 for pixels under thres (close to zero)
output = numerator./denominator;
output(tmp_mask) = 0;
end

