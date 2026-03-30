% Test for seeing how the fast code varies with slow code for the objective
% function either the Fredholm determinant or inverse of the Resolvent
% evaluated at two different points when we increase the number of points
% for a fixed geometry, i.e., no holes.

% chnkr.npt = nch*16;


zk1 = 3; zk2 = 4; 
rmax = 0.01:0.005:0.1; 
%ss = 1:0.05:5;
ss = 0.1;
NN = zeros(1,length(rmax));
rad = 1; ctr = [0.0;0.0];
circfun = @(t) ctr + rad*[cos(t(:).');sin(t(:).')];

nch = max(10, ceil(2*pi/rmax(1)) + 2);
optsc = [];
optsc.maxchunklen = ss;
prefc = []; 
prefc.nchmax = 10000;
chnkr = chunkerfuncuni(circfun, nch, optsc, prefc);
NN(1) = chnkr.npt;

tslow1 = zeros(1,length(rmax)); tfast1 = zeros(1,length(rmax));
tslow2 = zeros(1,length(rmax)); tfast2 = zeros(1,length(rmax));

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

for ii = 2:length(rmax) 

    nch = max(10, ceil(2*pi/rmax(ii)) + 2);
    opts = [];
    opts.maxchunklen = ss;
    pref = []; 
    pref.nchmax = 50000;
    chnkr = chunkerfuncuni(circfun, nch, opts, pref);    
    NN(ii) = chnkr.npt;

    if NN(ii)~=NN(ii-1)
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
save('fast_vs_slow_fixed_geometry.mat', 'tfast2',"tslow2","tfast1","tslow1", "NN");