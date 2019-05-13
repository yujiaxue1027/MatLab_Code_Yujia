% for i = 1:10
r = 2;
img = zeros(1024,1280,3);
img(512-r:512+r,640-r:640+r,:) = 255; 
fullscreen(uint8(img),2);
pause(10);
% end