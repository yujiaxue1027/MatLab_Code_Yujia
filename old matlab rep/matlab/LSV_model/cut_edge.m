function output = cut_edge(input,n)
output = input;
output(1:n,:,:) = 0;
output(end-(n-1):end,:,:) = 0;
output(:,1:n,:) = 0;
output(:,end-(n-1):end,:) = 0;
end

