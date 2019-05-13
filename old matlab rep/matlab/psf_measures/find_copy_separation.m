function [shift_length,shift_angle] = find_copy_separation(input)
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));
pad2d = @(x) padarray(x,0.5*size(x));
crop2d = @(x) x(1+size(x,1)/4:size(x,1)/4*3,...
    1+size(x,2)/4:size(x,2)/4*3);
reversed = @(x) flip(flip(x,1),2);
myxcorr2 = @(x,y) crop2d(Ft(F(pad2d(x)).*F(pad2d(reversed(y)))));
R = myxcorr2(input,input);
Rvec = R(:);
[~, indices] = sort(Rvec,'descend');
right_threshold = 500; % these two numbers are empirical
up_threshold = 100;
idx = 1;
find_it = 0;
[row0,col0] = ind2sub(size(R),indices(1));
while find_it ~= 1
    [row,col] = ind2sub(size(R),indices(idx));
    if (col-col0 > right_threshold) && (abs(row-row0) < up_threshold)
        right = col - col0;
        down = row - row0;
        find_it = 1;
    end
    idx = idx +1;
end
figure;
subplot(1,2,1), imagesc(input), axis image off;
subplot(1,2,2), imagesc(R), axis image off;
hold on;
plot(col0,row0,'r*');
plot(col0+right,row0+down,'r*');
shift_length = sqrt(right^2 + down^2);
shift_angle = rad2deg(atan(-down/right));
end

