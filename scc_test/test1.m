maxiter = 100;
rows = 256;
cols = 256;

for i = 1:maxiter
    filename = ['save/', num2str(i,'%.3d'),'.mat'];
    test_func(rows,cols,filename);
    disp(i);
end