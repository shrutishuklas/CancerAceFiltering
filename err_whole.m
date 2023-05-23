function [error_whole]=err_whole(y)
x=1:length(y);
sigma_xy = cumsum(x.*y);
sigma_x  = cumsum(x);
sigma_y  = cumsum(y);
sigma_xx = cumsum(x.*x);
n        = (1:length(y))';
det = n.*sigma_xx-sigma_x.*sigma_x;
mfwd = (n.*sigma_xy-sigma_x.*sigma_y)./det;
bfwd = -(sigma_x.*sigma_xy-sigma_xx.*sigma_y) ./det;

endpt=length(y);
error_whole= sum(abs((mfwd(endpt).*x(1:endpt)+bfwd(endpt))-y(1:endpt)));

end