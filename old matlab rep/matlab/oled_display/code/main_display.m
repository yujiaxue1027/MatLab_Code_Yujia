dataFolderName='inputImages/';
addpath(dataFolderName);
N=3; %images number
postpone=5; % postpone time / second
imgFormat='.jpg';
for i=1:N
    image2display=imread([dataFolderName,'image',num2str(i),imgFormat]);
    fullscreen(image2display,1);
    pause(postpone);
end
closescreen();
