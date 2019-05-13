%% compute NA of point source
data = [];
for i = 1:7
    img = im2double(imread([num2str(25+5*i),'.pgm']));
    data(:,:,i) = medfilt2(img,[5,5]);
end

i=7;
img = data(:,:,i);
max_val = max(max(img(401:1200,401:1200)));
[row,col] = find(img==max_val);
x = [1:2048]';
y = img(row,:)';
    
    