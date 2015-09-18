function [val,nspread]=nufft_not_good_kernel(x,oversamp,eps0)

if nargin<1, run_test; return; end;

nspread=floor(-log(eps0)/(pi*(oversamp-1)/(oversamp-.5)) + .5) + 1; 
nspread=nspread*2;
lambda=oversamp*oversamp * nspread/2 / (oversamp*(oversamp-.5));
tau=pi/lambda;
val=exp(-x.^2*tau).*(1-0.1*abs(x)/nspread); % a bit modified gaussian

end

function run_test

eps0=1e-4;
oversamp=2;
[~,nspread]=nufft_not_good_kernel(0,oversamp,eps0);
x=-nspread/2:0.01:nspread/2;
y=nufft_not_good_kernel(x,oversamp,eps0);
figure; plot(x,y);

optimize_nufft_kernel(@nufft_not_good_kernel,2);

end