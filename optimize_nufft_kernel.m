function optimize_nufft_kernel(kern,oversamp)

% This program tests a particular 1D kernel for use with nufft gridding.
% The first input (kern) is a function. See nufft_gaussian_kernel.m for an 
% example (also the default). The second input (oversamp) is the oversampling
% factor, which by default is 2. 
% 
% The function @kern will take three input parameters: (x, oversamp, eps0)
% x is the offset to compute the kern at. oversamp is the same as above, and
% eps0 is the desired precision. The output of @kern is [val,nspread] where 
% val is the value of the kernel and nspread is the width of the spreading
% kernel. Note that nspread is an OUTPUT parameter, not an INPUT parameter.
% In other words, @kern needs to do some optimization.
% 
% The program will generate two plots. Prescribed kernel size (nspread) as
% a function of prescribed epsilon, and the computed error as a function of 
% prescribed epsilon. Note that the user does not need to supply the
% correction function (ie, fourier transform of @kern) because it is computed
% automatically.
 
if nargin<1
    kern=@nufft_gaussian_kernel;
end

if nargin<2
    oversamp=2;
end;

N=360; % Theoretically it doesn't matter what we choose for N (I think)
epsilons=10.^(-(1:12)); % The range of epsilons to test

errs=[];
nspreads=[];
for ie=1:length(epsilons)
    eps0=epsilons(ie);
    [errs(ie),nspreads(ie)]=compute_max_error(N,oversamp,eps0,kern);
end;

figure; semilogx(epsilons,nspreads);
xlabel('Prescribed epsilon');
ylabel('Prescribed kernel size');

figure; loglog(epsilons,errs);
xlabel('Prescribed epsilon');
ylabel('Computed error');

end

function [err,nspread]=compute_max_error(N,oversamp,eps0,kern)

[~,nspread]=kern(0,oversamp,eps0);

NN=N*oversamp;
MM=ceil((NN+1)/2);
ns1=-ceil(nspread/2);
ns2=ns1+nspread-1;

A0=zeros(NN,1);
for j=ns1:ns2
    x0=MM;
    aa=round(x0);
    A0(aa+j)=kern((aa+j-x0),oversamp,eps0);
end;
B0=fftshift(fft(fftshift(A0)));

offsets=-1:0.01:1;
errs=zeros(1,length(offsets));
for io=1:length(offsets)
    offset=offsets(io);
    
    x0=MM+offset;
    aa=round(x0);
    
    A1=zeros(NN,1);
    for j=ns1:ns2
        A1(aa+j)=kern((aa+j-x0),oversamp,eps0);
    end;
    B1=fftshift(fft(fftshift(A1)));
    B1=B1./B0.*exp(-i*((1:NN)-MM)/NN*2*pi*offset)';
    %figure; plot(abs(B1-1));
    ind1=MM-ceil(N/2); ind2=ind1+N-1;
    errs(io)=max(abs(B1(ind1:ind2)-1));
end;

%figure; plot(offsets,errs);

err=max(errs);

end

