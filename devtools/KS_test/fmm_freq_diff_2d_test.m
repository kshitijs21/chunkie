clear;
%%%%%%%  Testing the "freq_diff" kernel at I don't know which level  %%%%%%%



%%%%%%%%   DEFINING THE CHUNKER FUNTION   %%%%%%%%%%%
% a = 2; b = 1;
% ge = @(t) [a*cos(t);b*sin(t)];
% der_ge = @(x) [-a*sin(x); b*cos(x)];
% der2_ge = @(x) [-a*cos(x); -b*sin(x)];
% nge = @(t) [b*cos(t);a*sin(t)]./( ((b*cos(t)).^2 + (a*sin(t)).^2).^0.5 );
% 
% 
% 
% jj = 100;
% pp = -pi : pi/jj : pi;
% dir = [cos(pp);sin(pp)];
% % dir = [cos(pp);sin(pp)];
% rng(1);                                    %%%%%% random number generator fixes seed
% qq = rand(1,length(pp));
% radr = 0.5 + 0.5*rand(1,length(pp));
% % sources = radr.*dir;
% sources = ((radr.').*(dir.')).';
% xp = sources(1,:); yp = sources(2,:);
% normalsrc = nge(pp);
% nxx = normalsrc(1,:); nyy = normalsrc(2,:);
% 
% %%%%%%%%%% CHARGES DEFINITION
% 
% charges1 = nxx.*qq;                         % normal_x.*charges
% charges2 = nyy.*qq;                         % normal_y.*charges 
% charges3 = xp.*charges1 + yp.*charges2;     % source_x.*normal_x.*charges + source_y.*normal_y.*charges
% 

%%%%%%%%%   SOURCES DEFINITION

% SOUCRE CORRESPONDING TO THE DOUBLE LAYER
% srcinfo = []; 
% srcinfo.sources = sources;
% srcinfo.dipstr = qq;
% srcinfo.dipvec = normalsrc;
% 
% srcinfo_kernel = [];
% srcinfo_kernel.r = sources;
% srcinfo_kernel.n = normalsrc;
% 
% % SOUCRE CORRESPONDING TO THE X-AXIS OF THE NORMAL FOR FREQ_DIFF
% srcinfo1 = []; 
% srcinfo1.sources = sources;
% srcinfo1.charges = charges1;
% 
% % SOUCRE CORRESPONDING TO THE Y-AXIS OF THE NORMAL FOR FREQ_DIFF
% srcinfo2 = []; 
% srcinfo2.sources = sources;
% srcinfo2.charges = charges2;
% 
% % SOUCRE CORRESPONDING TO THE COMPLETE BASIS FOR FREQ_DIFF
% srcinfo3 = []; 
% srcinfo3.sources = sources;
% srcinfo3.charges = charges3;
% 
% %%%%%%%%%%%% INITIALISING THE PARAMETERS
% 
% eps = 1e-11;
% pg = 1;
% zk = 1.41;
% h = 1e-2;
% hh = h./(1:1:100);
% deru1 = hfmm2d(eps, zk, srcinfo1, pg).pot;
% deru2 = hfmm2d(eps, zk, srcinfo2, pg).pot;
% deru3 = hfmm2d(eps, zk, srcinfo3, pg).pot;
% 
% deru = zk*xp.*deru1 + zk*yp.*deru2 - zk*deru3;
% deru2 = chnk.helm2d.fmm(eps, zk, srcinfo_kernel, srcinfo_kernel, 'freq_diff', qq);
% 
% for ii = 1:length(hh)
%     zp = zk + hh(ii);
%     zn = zk - hh(ii);
%     up = hfmm2d(eps, zp, srcinfo, pg).pot;
%     un = hfmm2d(eps, zn, srcinfo, pg).pot;
% 
%     err(ii) = norm((up - un)/hh(ii) - 2*deru2.');
%     err1(ii) = norm((up - un) - 2*hh(ii)*deru2.');
% end
% return

%%%%%%%  Testing the "freq_diff" kernel at layer potentials level  %%%%%%%

%%% Building the chunker  %%%

rad = 1; ctr = [0.0;0.0];
circfun = @(t) ctr + rad*[cos(t(:).');sin(t(:).')];

h = 4383;
rmax = 2/h;
nch = max(10, ceil(2*pi/rmax) + 2);  %%% number of panels

optss = [];
optss.maxchunklen = 1.1;
pref = []; 
pref.nchmax = 20000;

chnkr = chunkerfuncuni(circfun, nch, optss, pref);
% plot(chnkr, 'k.');

%%% Defining the sources %%%

% jj = chnkr.npt - 1;
% pp = 0 : pi/jj : pi;
% dir = [cos(pp);sin(pp)];
% rng(1);                      %%%%%% random number generator fixes seed
% qq = rand(1,length(pp));
% radr = 0.2 + 0.2*rand(1,length(pp));
% % sources = ((radr.').*(dir.')).';
% sources = dir;

% hold off
% plot(chnkr)
% hold on
% scatter(sources(1,:),sources(2,:),'x')
% % scatter(targets(1,:),targets(2,:),'x')
% axis equal 

%%% Obtaining the freq_diff %%%
xx = chnkr.r(1,:); yy = chnkr.r(2,:);
v = (xx - yy.^3 + xx.^2).*cos(-yy) + sin(yy).*(xx.^5 + xx);
zk = 2.41;

start1 = tic;
dkern = kernel('helmholtz','freq_diff',zk);
opts = [];
opts.corrections = true;
t1 = toc(start1);
% cormat = chunkermat(chnkr,dkern,opts);
cormat = chunkermat(chnkr, dkern, opts);
t1 = toc(start1);
start2 = tic;
opts = [];
opts.forcesmooth = true;
opts.cormat = cormat;
u_eval = chunkerkerneval(chnkr, dkern, v.', chnkr, opts);
t2 = toc(start2);
fprintf('nch = %d, t1 = %d, t2 = %d, t2/t1 = %d', nch, t1, t2, t2/t1)

return

h = 1e-3;
hh = h./(1:1:100);
for ii = 1:length(hh)
    zp = zk + hh(ii);
    zn = zk - hh(ii);
    %%%%%% Obtaining the vector v %%%
    Dp = kernel('helm', 'd', zp);  
    Ap = chunkermat(chnkr, Dp); 
    Dn = kernel('helm', 'd', zn);  
    An = chunkermat(chnkr, Dn);
    err(ii) = norm((Ap*(v.') - An*(v.'))/hh(ii) - 2*u_eval);
    %fprintf('i=%d,  cont=%d,  cent_diff=%d\n',ii, ...
    %       norm(Ap*(v.') - An*(v.')), norm((Ap*(v.') - An*(v.'))/hh(ii)));
    err1(ii) = norm(Ap*(v.') - An*(v.'));
end

% 
% figure, loglog(1:1:length(hh), err, 'k.');      % error is of order 100*h^2
% hold on; loglog(1:1:length(hh), hh.^2, 'r.'); 
% figure, loglog(1:1:length(hh), err1, 'm.');      % error is of order 100*h^2
% hold on; loglog(1:1:length(hh), hh, 'b.'); 
% figure, loglog(1:1:length(hh), err1, 'b.'); 
% hold on; loglog(1:1:length(hh), 100*hh.^3, 'm.');     % error is of order 100*h^3  

% Observations:
% 1. The rate of convergence in the Central Difference test is passed at
%    the Layer Potential level as well.
% 2. The time ratio of t_2/t_1 is of the order nch/4, where nch is the
%    number of panels on the chunker.