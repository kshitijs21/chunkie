% chunkermat_ostok2dTest0();

% function chunkermat_ostok2dTest0()

% CHUNKERMAT_OSTOK2DTEST
%
% test the matrix builder and do a basic solve

% % % % % iseed = 8675309;
% % % % % rng(iseed);
% % % % % 
% % % % % cparams = [];
% % % % % cparams.eps = 1.0e-10;
% % % % % cparams.nover = 1;
% % % % % pref = []; 
% % % % % pref.k = 16;
% % % % % narms = 3;
% % % % % amp = 0.25;
% start = tic;
% rad = 1; ctr = [0.0;0.0];
% circfun = @(t) ctr + rad*[cos(t(:).');sin(t(:).')];
% 
% 
% opts = [];
% opts.maxchunklen = 0.5;
% pref = []; pref.nchmax = 100000;
% chnkr1 = chunkerfunc(circfun, opts, pref);
% chnkr = chunkerfunc(circfun, opts, pref);
% 
% 
% % chnkr2 = 0.2*chnkr1;
% % chnkr2 = chnkr2.reverse;
% % chnkr = merge([chnkr1, chnkr2]); % only one hole
% 
% ns = 10;
% nt = 100;
% 
% sources = 3.05*circfun(1:ns);  % only one hole
% 
% 
% targets = 0.40*circfun(1:nt);
% 
% 
% strengths = randn(2*ns,1); % only one hole
% 
% 
% 
% sources_n = rand(2,ns);  % only one hole


iseed = 8675309;
% rng(iseed);

cparams = [];
cparams.eps = 1.0e-10;
cparams.nover = 1;
pref = []; 
pref.k = 20;
narms = 3;
amp = 0.25;
start = tic; chnkr = chunkerfunc(@(t) starfish(t,narms,amp),cparams,pref); 
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

% section %% 
zk = 0.3;
kernsp = kernel('ostok', 'sp', zk);
kerns = kernel('ostok', 's', zk);
kernd = kernel('ostok', 'd', zk);
kernsink = kernel('ostok', 'sinkv', zk);

% eval u on bdry


srcinfo = []; 
srcinfo.r = sources; 
srcinfo.n = sources_n;
kernmats = kerns.eval(srcinfo, chnkr); % for single layer
ubdry = kernmats*strengths; 

% kernmats = kernsink.eval(srcinfo, chnkr);
% kernmats = kernsp.eval(srcinfo, chnkr); % for sprime
% kernmats = kerns.eval(srcinfo, chnkr); % for single layer
% kernmats = kernd.eval(srcinfo, chnkr); % for double layer 
% ubdry = kernmats*strengths; 



% section %%
ubdry2 = kernd.eval(srcinfo, chnkr)*strengths;
n1 = chnkr1.npt;
ux1 = ubdry2(1:2:2*n1);
uy1 = ubdry2(2:2:2*n1);
nx1 = chnkr1.n(1,:).';
ny1 = chnkr1.n(2,:).';

r = sum((ux1.*nx1 + uy1.*ny1).*chnkr1.wts(:))


% section %%
% eval u at targets


targinfo = []; 
targinfo.r = targets;
targets_n = rand(2, nt); 
targets_n = targets_n./sqrt(targets_n(1,:).^2+targets_n(2,:).^2);
targinfo.n = targets_n;
kernmatstargs = kerns.eval(srcinfo, targinfo); % for sprime and single layer
% kernmatstargd = kernd.eval(srcinfo, targinfo); % for double layer 
utarg = kernmatstargs*strengths; 


%%%%% exact solution uisng analyticity not code
% section %%
% solve

% build Oscillatory Stokes dirichlet matrix
% fkern = kernel('ostok', 'sp', zk); % for sprime
fkernd = kernel('ostok', 'd', zk);  % for double layer
fkerns = kernel('ostok', 's', zk);  % for single layer
fkernsp = kernel('ostok', 'sp', zk);  % for single layer


% fkernd = kernel('ostok', 's', zk);  % for single layer
% uu = fkern.eval(chnkr, targinfo);


% section %%
start = tic; 
S = chunkermat(chnkr, fkerns);
D = chunkermat(chnkr, fkernd);
t1 = toc(start);

fprintf('%5.2e s : time to assemble matrix\n',t1)

% sys = 0.5*eye(size(D,1)) + D;
% sys = sys + normonesmat(chnkr)/sum(chnkr.wts(:)); % for sprime

sysd = -0.5*eye(size(D,1)) + D;
sysd = sysd + normonesmat(chnkr)/sum(chnkr.wts(:)); % for double layer

syssp = 0.5*eye(size(D,1)) + D;
syssp = syssp + normonesmat(chnkr)/sum(chnkr.wts(:)); % for double layer


% sys = D; % for single layer
syss = S ;

rhs = ubdry; 
rhs = rhs(:);

start = tic; 
sold = gmres(sysd,rhs,[],1e-12,1000); 
sols = gmres(syss,rhs,[],1e-12,1000); 
solsp = gmres(syssp,rhs,[],1e-12,1000); 
t1 = toc(start);

fprintf('%5.2e s : time for dense gmres\n',t1);

% evaluate at targets and compare

opts.usesmooth=false;
opts.verb=false;

%%%%%%%% SINGLE LAYER TEST 

fkerns =  kernel('ostok', 's', zk);
Ssol = chunkerkerneval(chnkr, fkerns, sols, targets, opts);
relerr = norm(utarg-Ssol,'fro')/(sqrt(chnkr.nch)*norm(utarg,'fro'));
fprintf('SL_ relative frobenius error %5.2e\n',relerr);

% fprintf('relative frobenius error %5.2e\n',relerr);
%return

%%%%%%%% DOUBLE LAYER TEST

fkernd =  kernel('ostok', 'd', zk);
Dsol = chunkerkerneval(chnkr, fkernd, sold, targets, opts);
relerr = norm(utarg-Dsol,'fro')/(sqrt(chnkr.nch)*norm(utarg,'fro'));
fprintf('DL_relative frobenius error %5.2e\n',relerr);

% fprintf('relative frobenius error %5.2e\n',relerr);

%%%%%%%   SPRIME TEST

fkernsp =  kernel('ostok', 'sp', zk); % for sprime
Spsol = chunkerkerneval(chnkr, fkernsp, solsp, targets, opts);
relerr = norm(utarg-Dsol,'fro')/(sqrt(chnkr.nch)*norm(utarg,'fro'));
fprintf('relative frobenius error %5.2e\n',relerr);
% 
% fprintf('relative frobenius error %5.2e\n',relerr);
% 
% assert(relerr < 1e-10);

%%%%%%%%%% SINK VELOCITY 

% end