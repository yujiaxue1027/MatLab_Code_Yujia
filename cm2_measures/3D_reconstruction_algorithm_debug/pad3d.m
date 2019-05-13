function output = pad3d(input)
    [r,c,~] = size(input);
    output = padarray(input,[r/2,c/2,0]);
end