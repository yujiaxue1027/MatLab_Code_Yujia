%% 888 is from grasshopper cam, 999 is lei's old data from ucb with PCO cam
%% 777 from ruilong's cam
img = imread('6.tif');
img = img(:);
img = sort(img);
img = diff(img);
img(img==0) = [];
min(img)