function [val,nspread,eps0]=nufft_gaussian_kernel(x,oversamp,param,opts)

if nargin<1, run_test; return; end;

eps0=10^(-(param-1));

nspread=floor(-2*log(eps0)/(pi*(oversamp-1)/(oversamp-.5)) + 3); 
lambda=oversamp*oversamp * nspread/2 / (oversamp*(oversamp-.5));
tau=pi/lambda;

val=exp(-x.^2*tau);


end

function run_test

param=5;
oversamp=2;

[~,nspread]=nufft_gaussian_kernel(0,oversamp,param);
x=-nspread/2:0.01:nspread/2;
y=nufft_gaussian_kernel(x,oversamp,param);
figure; plot(x,y);

end
