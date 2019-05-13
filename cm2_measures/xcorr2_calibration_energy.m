%% use autocorrelation to calibration energy of copies
figure,imagesc(gfp_xcorr),axis image;colormap(jet(256));colorbar
for i = 1:13
hold on;
plot(psf_xcorr_pos(i,2),psf_xcorr_pos(i,1),'r*');
end
xcorr_co = [];
for i = 1:13
    xcorr_co = [xcorr_co, gfp_xcorr(psf_xcorr_pos(i,1),psf_xcorr_pos(i,2))];
end 

fun_min_1 = @(x) (x(1)*x(7)+x(2)*x(8)+x(3)*x(9)-0.1298)^2 + ...
                 (x(1)*x(6)+x(4)*x(9)-0.1361)^2 + ...
                 (x(1)*x(5)+x(2)*x(6)+x(4)*x(8)+x(5)*x(9)-0.3901)^2 + ...
                 (x(1)*x(4)+x(2)*x(5)+x(3)*x(6)+x(4)*x(7)+x(5)*x(8)+x(6)*x(9)-0.5664)^2 + ...
                 (x(2)*x(4)+x(3)*x(5)+x(5)*x(7)+x(6)*x(8)-0.3906)^2 + ...
                 (x(3)*x(4)+x(6)*x(7)-0.1381)^2 + ...
                 (x(1)*x(3)+x(4)*x(6)+x(7)*x(9)-0.2302)^2 + ...
                 (x(1)*x(2)+x(2)*x(3)+x(4)*x(5)+x(5)*x(6)+x(7)*x(8)+x(8)*x(9)-0.6491)^2 + ...
                 (x(1)*x(9)-0.0237)^2 + ...
                 (x(1)*x(8)+x(2)*x(9)-0.0871)^2 + ...
                 (x(2)*x(7)+x(3)*x(8)-0.0898)^2 + ...
                 (x(3)*x(7)-0.0275)^2;

% (x(1)*x(1)+x(2)*x(2)+x(3)*x(3)+x(4)*x(4)+x(5)*x(5)+x(6)*x(6)+x(7)*x(7)+x(8)*x(8)+x(9)*x(9)-1)^2 + ...
             
             
x = fminsearch(fun_min_1,[0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3]);
    
    
    
    
    
    
    
    
    