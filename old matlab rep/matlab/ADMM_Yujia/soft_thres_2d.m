function [varargout] =  soft_thres_2d(v,h,tau)
    mag = sqrt(cat(1,v,zeros(1,size(v,2))).^2 + ...
            cat(2,h,zeros(size(h,1),1)).^2);
    magt = wthresh(mag,'s',tau);
    mmult = magt./mag;
    mmult(mag==0) = 0;
    varargout{1} = v.*mmult(1:end-1,:);
    varargout{2} = h.*mmult(:,1:end-1);
end
