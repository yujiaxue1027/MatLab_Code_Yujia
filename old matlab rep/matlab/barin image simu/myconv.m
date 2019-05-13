function output = myconv(img,kernel,method)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
channel = size(img,3);
output = zeros(size(img));
for i = 1:channel
    output(:,:,i) = conv2(img(:,:,i),kernel,method);
end





end

