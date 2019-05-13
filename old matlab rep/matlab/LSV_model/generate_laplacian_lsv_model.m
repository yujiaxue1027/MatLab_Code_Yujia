function PsiTPsi = generate_laplacian_lsv_model(rows,cols)

F = @(x) fftshift(fft2(ifftshift(x)));
% 
    lapl = zeros(rows,cols);    %Compute laplacian in closed form. This is the kernal to compute Psi'Psi
    rowpos = rows/2+1;
    colpos = cols/2+1;
    lapl(rowpos,colpos) = 4;
    
    lapl(rowpos,colpos-1) = -1;
    lapl(rowpos,colpos+1) = -1;    
    lapl(rowpos-1,colpos) = -1;
    lapl(rowpos+1,colpos) = -1;
    PsiTPsi = real(F(lapl));   %Compute power spectrum of laplacian
end