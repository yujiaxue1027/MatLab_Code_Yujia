function [output] = WT_operator(input)
average_hT = @(x) (x+circshift(x,[0,-1]))./sqrt(2);
difference_hT = @(x) (x-circshift(x,[0,-1]))./sqrt(2);
average_vT = @(x) (x+circshift(x,[-1,0]))./sqrt(2);
difference_vT = @(x) (x-circshift(x,[-1,0]))./sqrt(2);

% d: difference a: average v: vertical h:horizontal
dv = input(:,:,1);
dh = input(:,:,2);
av = input(:,:,3);
ah = input(:,:,4);
output = (difference_vT(dv)+difference_hT(dh)+...
    average_vT(av)+average_hT(ah))/4;
end
