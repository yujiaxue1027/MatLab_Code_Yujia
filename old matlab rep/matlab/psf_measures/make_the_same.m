function output = make_the_same(input,ref)
[refrow,refcol] = size(ref);
[inrow,incol] = size(input);
if inrow <= refrow
    output = zeros(refrow,refcol);
    row_offset = round((refrow - inrow)/2);
    col_offset = round((refcol - incol)/2);
    output(row_offset+1:row_offset+inrow,col_offset+1:col_offset+incol) = input;
else
    row_offset = round((inrow-refrow)/2);
    col_offset = round((incol-refcol)/2);
    output = input(row_offset+1:row_offset+refrow,col_offset+1:col_offset+refcol);
end
end

