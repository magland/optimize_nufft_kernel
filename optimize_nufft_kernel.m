%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Optimize the kernel for NUFFT Type 1 gridding in 1 dimension
function optimize_nufft_kernel

oversamp=2; % The oversampling factor can be changed
N=360; % Length of the uniform grid. Theoretically the choice of this parameter does not matter

nspreads=4:2:30; % Kernel spreading widths to test

% fac1s and fac2s are adjustment factors for optimizing the kernels
kb_fac1s=zeros(1,length(nspreads)); 
kb_fac2s=zeros(1,length(nspreads));
gaussian_fac1s=zeros(1,length(nspreads));

% These are the computed errors
kb_errs=zeros(1,length(nspreads));
gaussian_errs=zeros(1,length(nspreads));

FF1=figure; set(gcf,'Position',[150,50,600,600]);
FF2=figure; set(gcf,'Position',[800,50,600,600]);
for j=1:length(nspreads)
    nspread=nspreads(j);
    
    % Optimize the gaussian. Ideally val should be 1
    [VV,val]=fminsearch(@(V)objective_function_gaussian(V,nspread,oversamp,N),[1]);
    fprintf('gaussian: nspread: %d, fac1=%g, val=%g\n',nspread,VV(1),val);
    gaussian_fac1s(j)=VV(1);
    gaussian_errs(j)=10^val;
    
    % Optimize the Kaiser-Bessel kernel with two free parameters
    [VV,val]=fminsearch(@(V)objective_function_kb(V,nspread,oversamp,N),[1,1]);
    fprintf('kb: nspread: %d, fac1=%g, fac2=%g val=%g\n',nspread,VV(1),VV(2),val);
    kb_fac1s(j)=VV(1);
    kb_fac2s(j)=VV(2);
    kb_errs(j)=10^val;
    
    % Draw the results throughout
    figure(FF1);
    semilogy(nspreads,gaussian_errs,'b-o','LineWidth',3,'MarkerFaceColor','b','MarkerSize',6); hold on
    semilogy(nspreads,kb_errs,'r-o','LineWidth',3,'MarkerFaceColor','r','MarkerSize',6);
    xlabel('Kernel spread size'); ylabel('Computed precision');
    legend('Gaussian','Kaiser-Bessel');
    drawnow; hold off;
    
    % Draw the adjustment factors throughout
    figure(FF2);
    plot(nspreads,gaussian_fac1s,'b-o','LineWidth',3,'MarkerFaceColor','b','MarkerSize',6); hold on;
    plot(nspreads,kb_fac1s,'r-o','LineWidth',3,'MarkerFaceColor','r','MarkerSize',6);
    plot(nspreads,kb_fac2s,'m-o','LineWidth',3,'MarkerFaceColor','m','MarkerSize',6);
    xlabel('Kernel spread size'); ylabel('Optimization factor');
    legend('Gaussian factor 1','KB factor 1','KB factor 2');
    drawnow; hold off;
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Here is the Gaussian kernel
function [val,nspread,eps0]=nufft_gaussian_kernel(x,oversamp,nspread,opts)

fac1=opts.fac1;

%nspread=floor(-2*log(eps0)/(pi*(oversamp-1)/(oversamp-.5)) + 3); 
lambda=oversamp*oversamp * nspread/2 / (oversamp*(oversamp-.5));
tau=pi/lambda*fac1;

val=exp(-x.^2*tau);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Here is the Kaisser-Bessel kernel
function [val,nspread]=nufft_kb_kernel(x,oversamp,nspread,opts)

W=nspread*opts.fac1;
beta=pi*sqrt(W*W/4-0.8)*opts.fac2;

y=beta*sqrt(1-(2*x/W).^2);
val=besseli(0,y);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The objection for KB, to be minimized
function ret=objective_function_kb(V,nspread,oversamp,N)

kern=@nufft_kb_kernel;
kern_opts.fac1=V(1);
kern_opts.fac2=V(2);

nspreads=nspread;
err=optimize_nufft_kernel_2(kern,kern_opts,nspreads,oversamp,N);
ret=log(err)/log(10);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The objective function for Gaussian, to be minimized
function ret=objective_function_gaussian(V,nspread,oversamp,N)

kern=@nufft_gaussian_kernel;
kern_opts.fac1=V(1);

nspreads=nspread;
err=optimize_nufft_kernel_2(kern,kern_opts,nspreads,oversamp,N);
ret=log(err)/log(10);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part of the objective function
function [errs,nspreads]=optimize_nufft_kernel_2(kern,kern_opts,nspreads,oversamp,N)

errs=[];
for ip=1:length(nspreads)
    nspread0=nspreads(ip);
    [errs(ip)]=compute_max_error(N,oversamp,nspread0,kern,kern_opts);
end;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Here we compute the maximum error numerically
function [err,nspread]=compute_max_error(N,oversamp,nspread0,kern,kern_opts)

[~,nspread]=kern(0,oversamp,nspread0,kern_opts);

NN=N*oversamp;
MM=ceil((NN+1)/2);
ns1=-ceil(nspread/2);
ns2=ns1+nspread-1;

A0=zeros(NN,1);
for j=ns1:ns2
    x0=MM;
    aa=round(x0);
    A0(aa+j)=kern((aa+j-x0),oversamp,nspread0,kern_opts);
end;
B0=fftshift(fft(fftshift(A0))); % This the calibration term - for correction

offsets=-0.5:0.05:0.5; % Try different offsets
errs=zeros(1,length(offsets));
for io=1:length(offsets)
    offset=offsets(io);
    
    x0=MM+offset;
    aa=round(x0);
    
    A1=zeros(NN,1);
    for j=ns1:ns2
        A1(aa+j)=kern((aa+j-x0),oversamp,nspread0,kern_opts);
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

