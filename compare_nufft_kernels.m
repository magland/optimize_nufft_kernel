addpath([fileparts(mfilename('fullpath')),'/colorspace']);

list={};

fac2=1.6;

kb_fac1=[0.75, 0.80, 0.85, 0.90, 0.95, 1.00];
kb_fac2=[fac2, fac2, fac2, fac2, fac2, fac2];

colors=distinguishable_colors(30,{'k','w'});

for jj=1:length(kb_fac1)
    A=struct;
    A.kern=@nufft_trivial_kernel;
    A.opts.fac1=kb_fac1(jj); A.opts.fac2=kb_fac2(jj);
    A.name=sprintf('Kaiser Bessel fac1=%g, fac2=%g',kb_fac1(jj),kb_fac2(jj));
    A.color=colors(1+mod(jj-1,length(colors)),:);
    list{end+1}=A;
end;

A=struct;
A.kern=@nufft_gaussian_kernel;
A.opts=struct;
A.name='Gaussian';
A.color='k';
list{end+1}=A;

params=1:10;
oversamp=3;
N=360;

figure; set(gcf,'position',[0,0,800,800]);
legends={};
for ii=1:length(list)
    kern=list{ii}.kern;
    opts=list{ii}.opts;
    name=list{ii}.name;
    [errs,nspreads]=optimize_nufft_kernel(kern,opts,params,oversamp,N);
    hh=semilogy(nspreads,errs,'LineWidth',3);
    set(hh,'Color',list{ii}.color);
    legends{end+1}=name;
    hold on;
end;
xlabel('Kernel spread size');
ylabel('Computed precision');
legend(legends);


%legend('Gaussian','Kaiser-Bessel');
