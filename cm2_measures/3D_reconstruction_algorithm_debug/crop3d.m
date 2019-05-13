function output = crop3d(input)
    [r,c,z] = size(input);
    output = input(r/4+1:r*3/4,c/4+1:c*3/4,z/2);
end