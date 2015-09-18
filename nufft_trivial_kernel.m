function [val,nspread]=nufft_trivial_kernel(x,oversamp,eps0)

nspread=10;

W=nspread;
beta=pi*sqrt(W*W/4-0.8);
y=beta*sqrt(1-(2*x/W).^2);
val=besseli(0,y);

end

