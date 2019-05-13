function output_pattern = generate_calibration_pattern(square_size,copies)
%% generate calibration pattern for LSV PSF calibration
% square_size = 32;
radius = square_size/2/sqrt(3)*2;
[x,y] = meshgrid(-square_size:1:square_size-1);
pattern_large = 255*ones(square_size);
pattern_large = padarray(pattern_large,0.5*size(pattern_large));
c1 = [-1*(square_size/2-radius/2),0]';
for i = square_size/2+1 : square_size/2 + square_size
    for j = 1:square_size/2
        pos = [x(i,j),y(i,j)]';
        dist = sqrt(sum((c1-pos).^2));
        if dist <= radius
            pattern_large(i,j) = 255;
        end
    end
end
c2 = [+1*(square_size/2-radius/2),0]';
for i = square_size/2+1 : square_size/2 + square_size
    for j = square_size/2 + square_size:square_size*2
        pos = [x(i,j),y(i,j)]';
        dist = sqrt(sum((c2-pos).^2));
        if dist <= radius
            pattern_large(i,j) = 255;
        end
    end
end
c3 = [0,-1*(square_size/2+radius/2)]';
for i = 1+square_size/2:square_size
    for j = square_size/2+1 : square_size/2 + square_size
        pos = [x(i,j),y(i,j)]';
        dist = sqrt(sum((c3-pos).^2));
        if dist <= radius
            pattern_large(i,j) = 0;
        end
    end
end
c4 = [0,+1*(square_size/2+radius/2)]';
for i = 1+square_size:square_size/2 + square_size
    for j = square_size/2+1 : square_size/2 + square_size
        pos = [x(i,j),y(i,j)]';
        dist = sqrt(sum((c4-pos).^2));
        if dist <= radius
            pattern_large(i,j) = 0;
        end
    end
end
% figure,imagesc(pattern_large),axis image off,colormap gray;truesize;

%%
% copies = 5;
basis = pattern_large;
outsize = copies*size(basis);
output = zeros(outsize);
pattern = padarray(basis,0.5*(size(output)-size(basis)));
[x,y] = meshgrid(-floor(copies/2):1:-floor(copies/2)+copies-1);
shift_basis1 = [0,2*square_size];
shift_basis2 = [square_size,square_size];
for i = 1:size(x,1)
    for j = 1:size(x,2)
        shift = x(i,j)*shift_basis1 + y(i,j)*shift_basis2;
        output = output + circshift(pattern,shift);
    end
end
valid_size = copies*square_size;
output = output(size(output,1)/2-valid_size/2+1:size(output,1)/2+valid_size/2,...
    size(output,2)/2-valid_size/2+1:size(output,2)/2+valid_size/2);
% figure,imagesc(output),axis image off,colormap gray;truesize;
output_pattern = output;