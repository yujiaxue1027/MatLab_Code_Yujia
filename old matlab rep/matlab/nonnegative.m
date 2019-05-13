function [ output ] = nonnegative( input )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
output = input;
idx = input<0;
output(idx) = 0;

end

