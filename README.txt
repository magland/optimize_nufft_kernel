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