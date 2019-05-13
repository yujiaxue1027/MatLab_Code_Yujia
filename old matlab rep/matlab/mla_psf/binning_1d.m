function output = binning_1d(input, binning)
input = input(:);
num_in = length(input);
num_out = ceil(num_in/binning);
output = zeros(num_out,1);
num_in_new = num_out*binning;
input = padarray(input,[num_in_new-num_in,0],'post');
for i = 1:num_out
    output(i) = mean(input(1+(i-1)*binning:i*binning));
end

