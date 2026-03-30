% Test for seeing how the fast code varies with slow code for the objective
% function either the Fredholm determinant or inverse of the Resolvent
% evaluated at two different points when we increase the number of holes
% by changing the geometry, i.e., introducing more number of holes.

zk1 = 3; zk2 = 4; 
xx = 1:0.1:2;
nh = zeros(1,length(xx));   % number of holes for get_chunker(xx(ii))
NN = zeros(1,length(xx));   % chnkr.npt for get_chunker(xx(ii))
[chnkr, nh(1)] = get_chunker(xx(1));
NN(1) = chnkr.npt;
tslow1 = zeros(1,length(xx)); tfast1 = zeros(1,length(xx));
tslow2 = zeros(1,length(xx)); tfast2 = zeros(1,length(xx));

func = 'fred';

opts = [];
opts.speed = 'slow';

t1 = tic;
f1s = fred_opti(zk1, chnkr, func, opts);
tslow1(1,1) = toc(t1);

opts1 = [];
opts1.speed = 'fast';

t2 = tic;
f1f = fred_opti(zk1, chnkr, func, opts1);
tfast1(1,1) = toc(t2);


t1 = tic;
f2s = fred_opti(zk2, chnkr, func, opts);
tslow2(1,1) = toc(t1);


t2 = tic;
f2f = fred_opti(zk2, chnkr, func, opts1);
tfast2(1,1) = toc(t2);

for ii = 2:length(xx)    
    [chnkr,nh(ii)] = get_chunker(xx(ii));
    NN(ii) = chnkr.npt;

    if nh(ii)~=nh(ii-1)
        t1 = tic;
        f1s = fred_opti(zk1, chnkr, func, opts);
        tslow1(ii) = toc(t1);
        t2 = tic;
        f1f = fred_opti(zk1, chnkr, func, opts1);
        tfast1(ii) = toc(t2);

        t11 = tic;
        f2s = fred_opti(zk2, chnkr, func, opts);
        tslow2(ii) = toc(t11);
        t21 = tic;
        f2f = fred_opti(zk2, chnkr, func, opts1);
        tfast2(ii) = toc(t21);
    else
      
        tslow1(ii) = tslow1(ii-1);
        tfast1(ii) = tfast1(ii-1);
        tslow2(ii) = tslow2(ii-1);
        tfast2(ii) = tfast2(ii-1);
    end
end
save('fast_vs_slow_varying_geometry_l_from_1_to_2.mat', 'tfast2',"tslow2","tfast1","tslow1", "NN");