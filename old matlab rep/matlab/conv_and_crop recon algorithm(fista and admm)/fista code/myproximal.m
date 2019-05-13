function [ output ] = myproximal( input,thres)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[rows,cols] = size(input);
input = double(input);
% gx = input - [zeros(rows,1),input(:,1:end-1)];
% gy = input - [zeros(1,cols);input(1:end-1,:)];
% apply soft thresholding in gradient magnitude, does not work due to
% dividing by zero
% gmag = sqrt(gx.^2+gy.^2);
% gratio = (gx+eps)./(gy+eps);
% new_gmag = wthresh(gmag,'s',thres);
% new_gx = new_gmag.*gratio./(sqrt(1+gratio.^2));
% new_gy = new_gmag./(sqrt(1+gratio.^2));

% apply soft thresholding in each direction with thres = thres./sqrt of 2
% does not work
% new_gx = wthresh(gx,'s',thres/sqrt(2));
% new_gy = wthresh(gy,'s',thres/sqrt(2));
% rex = cumsum(new_gx,2);
% rey = cumsum(new_gy,1);
% output = (rex+rey)./2;

[x,y] = wavedec2(input,1,'haar');
x = wthresh(x,'s',thres);
output = waverec2(x,y,'haar');
end

