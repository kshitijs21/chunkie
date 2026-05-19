clear all;

iseed = 8675309;
rng(iseed);

cparams = [];
cparams.eps = 1.0e-10;
cparams.nover = 1;
pref = []; 
pref.k = 20;
narms = 3;
amp = 0.25;
start = tic; 
chnkr = chunkerfunc(@(t) starfish(t,narms,amp),cparams,pref); 
t1 = toc(start);

fprintf('%5.2e s : time to build geo\n',t1);

% sources

ns = 10;
ts = 0.0+2*pi*rand(ns,1);
sources = starfish(ts,narms,amp);
sources = 3.0*sources;
strengths = randn(2*ns,1);
sources_n = rand(2,ns);

% targets

nt = 100;
ts = 0.0+2*pi*rand(nt,1);
targets = starfish(ts,narms,amp);
targets = targets.*repmat(rand(1,nt),2,1)*0.8;

plot(chnkr, 'r.'); hold on;
plot(targets(1,:), targets(2,:), 'kx')
hold on;
plot(sources(1,:), sources(2,:), 'bo')
axis equal




targs = chnkr.r; targs = reshape(targs,2,chnkr.k*chnkr.nch);
targstau = tangents(chnkr); 
targstau = reshape(targstau,2,chnkr.k*chnkr.nch);

plot(chnkr, 'r.'); hold on;
plot(targets(1,:), targets(2,:), 'kx');
plot(sources(1,:), sources(2,:), 'bo');hold on;

t1 = toc(start);

fprintf('%5.2e s : time to build geo\n',t1)


zk = 0.3;


% eval u on bdry


srcinfo = []; 
srcinfo.r = sources; 
srcinfo.n = sources_n;

targinfo = []; 
targinfo.r = targets;
targets_n = rand(2, nt); 
targets_n = targets_n./sqrt(targets_n(1,:).^2+targets_n(2,:).^2);
targinfo.n = targets_n;


srch = srcinfo.r(:,:);
targh = targinfo.r(:,:);
[val,grad,hess,der3, ~,g0,g1,g21,g3,g4] = chnk.obihar2d.green(zk,srch,targh);
c = g0;
[valf,gradf,hessf,der3f,~,~,gf0,gf1,gf21,gf3,gf4] = chnk.flex2d.helmdiffgreen(zk,srch,targh);
cc = gf0;
valf = valf/(zk^2);
gradf = gradf/(zk^2);
hessf = hessf/(zk^2);
der3f = der3f/(zk^2);

fprintf('KS_test_val = %d\n', norm(val(:)-valf(:)));

fprintf('KS_test_grad1 = %d\n', norm(grad(:,:,1)-gradf(:,:,1)));
fprintf('KS_test_grad2 = %d\n', norm(grad(:,:,2)-gradf(:,:,2)));

fprintf('KS_test_hess1 = %d\n', norm(hess(:,:,1)-hessf(:,:,1)));
fprintf('KS_test_hess2 = %d\n', norm(hess(:,:,2)-hessf(:,:,2)));
fprintf('KS_test_hess3 = %d\n', norm(hess(:,:,3)-hessf(:,:,3)));

fprintf('KS_test1_der3 = %d\n', norm(der3(:,:,1)-der3f(:,:,1)));
fprintf('KS_test2_der3 = %d\n', norm(der3(:,:,2)-der3f(:,:,2)));
fprintf('KS_test3_der3 = %d\n', norm(der3(:,:,3)-der3f(:,:,3)));
fprintf('KS_test3_der4 = %d\n', norm(der3(:,:,4)-der3f(:,:,4)));


fprintf('KS_test1_g0 = %d\n', norm(g0(:)-gf0(:)));
fprintf('KS_test2_g1 = %d\n', norm(g1(:)-gf1(:)));
fprintf('KS_test3_g21 = %d\n', norm(g21(:)-gf21(:)));
fprintf('KS_test3_g3 = %d\n', norm(g3(:)-gf3(:)));
fprintf('KS_test3_g4 = %d\n', norm(g4(:)-gf4(:)));
