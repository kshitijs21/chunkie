%l = [1.5 2 2.6];
l = [3];
z0 = 3;

func = 'fred';
method = 'cheb';
%kk = 0.001:0.0045:0.01;
kk = 0.01;
cheb1 = z0-kk;
cheb2 = z0+kk;
opts = [];

tt = zeros(1,length(kk));
deg = zeros(1,length(kk));

for ii = 1:length(l)
    [chnkr,nh] = get_chunker(l(ii));
    ff = @(zk) fred_opti(zk, chnkr, func, opts);
    fprintf('Hole = %d:\n', nh);
    for jj = 1:length(kk)
        
        a = cheb1(jj);
        b = cheb2(jj);
        t1 = tic;
        f = chebfun(ff, [a, b]);
        tt(jj) = toc(t1);
        deg(jj) = length(f);
        fprintf('  Interval = [%d, %d], Time by Chebfun = %d, Degree of Interpolant = %d, \n', a, b,  tt(jj), deg(jj));
    end
    save(['degree_vs_hole_plot_for_' num2str(nh) '_hole.mat'], 'cheb1','cheb2', 'tt', 'deg');
end