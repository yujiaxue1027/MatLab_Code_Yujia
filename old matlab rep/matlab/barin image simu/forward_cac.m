function output = forward_cac(x,psf,length)
output = myconv(x,psf,'same');
tmp_size = size(output,1);
output = output(round(tmp_size/2)-length/2:round(tmp_size/2)+length/2-1,...
    round(tmp_size/2)-length/2:round(tmp_size/2)+length/2-1);
end

