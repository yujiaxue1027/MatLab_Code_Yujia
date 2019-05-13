function [output] = W_operator(img)
average_h = @(x) (x+circshift(x,[0,1]))./sqrt(2);
difference_h = @(x) (x-circshift(x,[0,1]))./sqrt(2);
average_v = @(x) (x+circshift(x,[1,0]))./sqrt(2);
difference_v = @(x) (x-circshift(x,[1,0]))./sqrt(2);


% d: difference a: average v: vertical h:horizontal
dv = difference_v(img);
dh = difference_h(img);
av = average_v(img);
ah = average_h(img);
output = [];
output(:,:,1) = dv;
output(:,:,2) = dh;
output(:,:,3) = av;
output(:,:,4) = ah;

end
