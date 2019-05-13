function [output] = tv_proximal(input,thres)
    input = double(input);
    [x,y] = wavedec2(input,1,'haar');
    x = wthresh(x,'s',thres);
    output = waverec2(x,y,'haar');
end