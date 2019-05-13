function output = softthres_d(input,thres)
output = input;
output(:,:,1:2) = wthresh(output(:,:,1:2),'s',thres);
end