rad = 1; ctr = [0.0;0.0];
circfun = @(t) ctr + rad*[cos(t(:).');sin(t(:).')];

nch = 500;
% nch = max(10, ceil(2*pi/rmax) + 2);
optsc = [];
optsc.maxchunklen = 0.02;
prefc = []; 
prefc.nchmax = 10000;
chnkr = chunkerfuncuni(circfun, nch, optsc, prefc);

Dk = kernel('helm', 'd', 3); tic, A=chunkermat(chnkr, Dk); toc;
tic, cc = det(A); toc