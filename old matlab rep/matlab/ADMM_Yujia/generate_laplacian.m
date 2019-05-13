function PsiTPsi = generate_laplacian(rows,cols)

F = @(x) fftshift(fft2(ifftshift(x)));
% 
    lapl = zeros(2*rows,2*cols);    %Compute laplacian in closed form. This is the kernal to compute Psi'Psi
    lapl(rows+1,cols+1) = 4;
    
    lapl(rows+1,cols+1+1) = -1;
    lapl(rows+1+1,cols+1) = -1;
    lapl(rows,cols+1) = -1;
    lapl(rows+1,cols) = -1;
    PsiTPsi = abs(F(lapl));   %Compute power spectrum of laplacian
end