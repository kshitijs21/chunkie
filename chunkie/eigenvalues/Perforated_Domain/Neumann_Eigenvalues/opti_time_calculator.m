[chnkr,nh] = get_chunker(5);
xx = chnkr.r(1,:);
yy = chnkr.r(1,:);
uu = @(x,y) (21 + 4*x + 635*y + 6782*y.^2).*sin(354 + x + 23*x.*y);
vv = @(x,y) (21 + 635*x.*y + 6782*y.^4).*cos(354 + y.^2 + 23*x.*(y.^2));
u = uu(xx,yy);
v = vv(xx,yy);
zk = 3;

%%%%%%%%    slow code    %%%%%%%%%%%

t0 = tic;
Dk = 2*kernel('helm', 'd', zk); 
timeek = toc(t0);
t1 = tic;
A = chunkermat(chnkr, Dk);    
timee = toc(t1);
t2 = tic;
A = A + eye(chnkr.npt);
timee1 = toc(t2);
t3 = tic;
opti = sum( ones(1, length(u))*( ( (u.').*( A\(v.') ) ).*(chnkr.wts(:)) ) );
timee2 = toc(t3);
t4 = tic;
f = 1/opti;
timee3 = toc(t4);
fprintf('Hole = %d, \n kernel buildup time = %d,\n chnukermat = %d,\n id sum = %d,\n opti = %d,\n f = %d\n', nh, timeek, timee, timee1, timee2, timee3);

%%%%%%%%    fast code    %%%%%%%%%%%

opts_flam = [];
opts_flam.flamtype = 'rskelf';
opts_flam.forceproxy = true;
opts_flam.occ = 200;


t01 = tic;
F = chunkerflam(chnkr, Dk, 1.0, opts_flam);
timeek1 = toc(t01);
t11 = tic;
sol = (rskelf_sv(F,v.'));
timee01 = toc(t11);
t21 = tic;
opp = (u.').*sol.*chnkr.wts(:);
timee11 = toc(t21);
t31 = tic;
opti = sum(opp(:));  
timee21 = toc(t31);
t41 = tic;
f = 1/opti;
timee31 = toc(t41);
fprintf('Hole = %d, \n kernel buildup time = %d,\n chnukermat = %d,\n id sum = %d,\n opti = %d,\n f = %d\n', nh, timeek1, timee01, timee11, timee21, timee31);