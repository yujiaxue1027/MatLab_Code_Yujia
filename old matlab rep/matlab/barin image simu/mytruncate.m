function [ output ] = mytruncate( input ,vmin,vmax)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
output = input;
idx = input<vmin;
output(idx) = vmin;
idx = input>vmax;
output(idx) = vmax;

end

