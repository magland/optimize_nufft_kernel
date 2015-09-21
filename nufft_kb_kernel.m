function [val,nspread]=nufft_trivial_kernel(x,oversamp,param0,opts)

nspread=param0*2;

W=nspread*opts.fac1;
beta=pi*sqrt(W*W/4-0.8)*opts.fac2;

y=beta*sqrt(1-(2*x/W).^2);
val=besseli(0,y);

end

