function pcc = compute_pcc(img1,img2)
img1 = img1(:) - mean(img1(:));
img2 = img2(:) - mean(img2(:));
numerator = mean(img1.*img2);
denominator = sqrt(mean(img1.^2)*mean(img2.^2));
pcc = numerator/denominator;
end

