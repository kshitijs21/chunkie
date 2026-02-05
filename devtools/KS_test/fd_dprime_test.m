clear;

%%%%%%%%%% TEST FOR fd_dprime kernel %%%%%%%%%%%

chnkr = get_chunker(1);
hh = 1:1:40;
err = zeros(1,length(hh));
h0 = 1e-2;
zz = 2.8882;
Freqk = kernel('helm', 'fd_dprime', zz);
Ak = chunkermat(chnkr, Freqk); 
for i = 1:length(hh)
     h = h0/hh(i);
     Dp = kernel('helm', 'dprime', zz + h);  
     Ap = chunkermat(chnkr, Dp);
     Dm = kernel('helm', 'dprime', zz - h);
     Am = chunkermat(chnkr, Dm);
     err(i) = norm((1/h)*(Ap - Am) - 2*Ak);
end
loglog(hh, err, 'm.');
hold on; loglog(hh, (h0./hh).^2, 'b.');

%%%%%%%%%%  FMM test  %%%%%%%%%%%%%%

rad = 1; ctr = [0.0;0.0];
circfun = @(t) ctr + rad*[cos(t(:).');sin(t(:).')];

hl = 4;
rmax = 2/hl;
nch = max(10, ceil(2*pi/rmax) + 2);  %%% number of panels

optss = [];
optss.maxchunklen = 1.1;
pref = []; 
pref.nchmax = 20000;

chnkrl = chunkerfuncuni(circfun, nch, optss, pref);

%%% Obtaining the fd_dprime %%%
xx = chnkrl.r(1,:); yy = chnkrl.r(2,:);
v = (xx - yy.^3 + xx.^2).*cos(-yy) + sin(yy).*(xx.^5 + xx);
zk = 2.41;

start1 = tic;
dkernl = kernel('helmholtz','fd_dprime',zk);
% Fkern = chunkermat(chnkr, dkern); 
opts = [];
opts.corrections = true;
t1 = toc(start1);
% cormat = chunkermat(chnkr,dkern,opts);
cormat = chunkermat(chnkrl, dkernl, opts);
t1 = toc(start1);
fprintf('Time to compute the weights=%d\n', t1);
start2 = tic;
opts = [];
opts.forcesmooth = true;
opts.cormat = cormat;
u_eval = chunkerkerneval(chnkrl, dkernl, v.', chnkrl, opts);
t2 = toc(start2);
fprintf('Time to evaluate the kernel at the boundary level=%d\n', t2);
% return

hl = 1e-2;
hhl = hl./(1:1:100);
errl = zeros(1,length(hhl));
for ii = 1:length(hhl)
    zpl = zk + hhl(ii);
    znl = zk - hhl(ii);
    Dpl = kernel('helm', 'dprime', zpl);  
    Apl = chunkermat(chnkrl, Dpl); 
    Dnl = kernel('helm', 'dprime', znl);  
    Anl = chunkermat(chnkrl, Dnl);
    errl(ii) = norm((Apl*(v.') - Anl*(v.'))/hhl(ii) - 2*u_eval);
end

figure, loglog(1:1:length(hhl), errl, 'k.');      
hold on; loglog(1:1:length(hhl), hhl.^2, 'r.'); 