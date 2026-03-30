clear;
% Interpolating the optimization function requires around 300 degree
% polynomials
% ncheb = 30;

% for ii = 1:10
%     a = rand(1000*ii);
%     t1 = tic;
%     dd = det(a);
%     ttt(ii) = toc(t1);
% end
% loglog(1:1:10, ttt, 'k.');
% return

% rad = 1; ctr = [0.0;0.0];
% circfun = @(t) ctr + rad*[cos(t(:).');sin(t(:).')];
% ss = 0.1:0.05:5;
% NN = zeros(1,length(ss));
% l = 3; k = 1/(4*l); h = 5*l; 
% rmax = 0.01:0.005:1;
% 
% for ii = 1:length(ss)
%     nch = max(10, ceil(2*pi/rmax(ii)) + 2);
%     opts = [];
%     opts.maxchunklen = ss(ii);
%     pref = []; 
%     pref.nchmax = 10000*ii;
% 
%     chnkr = chunkerfuncuni(circfun, nch, opts, pref);
%     NN(ii) = chnkr.npt;
% end
% plot(NN, 'k.');
% 
% return


% chnkr = get_chunker(1);
% xx = chnkr.r(1,:);
% yy = chnkr.r(2,:);
% u = (xx.^2 - yy.^2 + yy).*sin(xx) + cos(xx).*(0.5*yy.^2 + 0.3*yy);
% v = (xx - yy.^3 + xx.^2).*cos(-yy) + sin(yy).*(xx.^2 + xx);
% A = @(zk) get_matrix(zk, chnkr);
% g = @(zk) ones(1, length(u))*( ( (u.').*( A(zk)\(v.') ) ).*(chnkr.wts(:)) );
% a = 3.3; b = 3.5; deg =  25;
% ff = chebfun(@(zk) 1/g(zk), [a,b], deg);

% opts = [];
% opts.u = u;
% opts.v = v;
% tic, fred_opti(zk, chnkr, 'ir', opts); toc

zk1 = 3; zk2 = 4; 
% rmax = 0.01:0.005:1; 
% ss = 1:0.05:5;
% NN = zeros(1,length(ss));
rmax = 0.25;

rad = 1; ctr = [0.0;0.0];
circfun = @(t) ctr + rad*[cos(t(:).');sin(t(:).')];

nch = 500;
% nch = max(10, ceil(2*pi/rmax) + 2);
optsc = [];
optsc.maxchunklen = 0.02;
prefc = []; 
prefc.nchmax = 10000;
chnkr = chunkerfuncuni(circfun, nch, optsc, prefc);
NN = chnkr.npt;


% a = 3; b = 3.645; deg =  17;
% xx = 1:0.1:4;
% nh = zeros(1,length(xx)); 
%N = zeros(1,length(xx));
% tslow1 = zeros(1,length(ss)); tfast1 = zeros(1,length(ss));
% tslow2 = zeros(1,length(ss)); tfast2 = zeros(1,length(ss));


opts = [];
opts.speed = 'slow';

t1 = tic;
f1s = fred_opti(zk1, chnkr, 'fred', opts);
tslow1(1,1) = toc(t1);

opts1 = [];
opts1.speed = 'fast';

t2 = tic;
f1f = fred_opti(zk1, chnkr, 'fred', opts1);
tfast1(1,1) = toc(t2);


t1 = tic;
f2s = fred_opti(zk2, chnkr, 'fred', opts);
tslow2(1,1) = toc(t1);


t2 = tic;
f2f = fred_opti(zk2, chnkr, 'fred', opts1);
tfast2(1,1) = toc(t2);
% N(1) = chnkr.npt;

return

for ii = 2:length(ss)

    nch = max(10, ceil(2*pi/rmax(ii)) + 2);
    opts = [];
    opts.maxchunklen = ss(ii);
    pref = []; 
    pref.nchmax = 10000*ii;
    chnkr = chunkerfuncuni(circfun, nch, opts, pref);
    
    % [chnkr,nh(ii)] = get_chunker(xx(ii));
    NN(ii) = chnkr.npt;
    if NN(ii)~=NN(ii-1)
        t1 = tic;
        f1s = fred_opti(zk1, chnkr, 'fred', opts);
        tslow1(ii) = toc(t1);
        t2 = tic;
        f1f = fred_opti(zk1, chnkr, 'fred', opts1);
        tfast1(ii) = toc(t2);

        t11 = tic;
        f2s = fred_opti(zk2, chnkr, 'fred', opts);
        tslow2(ii) = toc(t11);
        t21 = tic;
        f2f = fred_opti(zk2, chnkr, 'fred', opts1);
        tfast2(ii) = toc(t21);
    else
      
        tslow1(ii) = tslow1(ii-1);
        tfast1(ii) = tfast1(ii-1);
        tslow2(ii) = tslow2(ii-1);
        tfast2(ii) = tfast2(ii-1);
    end
end
return

%%
% plot(imag(ff), 'k.');
% rts = roots(ff, 'complex');
% 
% 

%%
%Amat = get_matrix(akk, chnkr); 

function [A] = get_matrix(zk, chnkr1)

    Dk = 2*kernel('helm', 'd', zk);  
    A = chunkermat(chnkr1, Dk);   
%%
    A = A + eye(chnkr1.npt);
end


%%%%%%%%% OBSERVATIONS FROM THIS CODE %%%%%%%%
% 1. For the number of holes = 12 and for the function 
%    ff = chebfun(@(zk) 1/g(zk), [3.4,3.5]); 
%    we see that only 18 degree polynomials is needed,
%    which means that if we have made precise the interval where 
%    we are looking for the eigenvalues then  we are able to use
%    this function as the optimization function and try it
%    to calculate the eigenvalue of the operator.
%    Also, ur agrees with the Fredholm determinant technique
% 2. For the number of holes = 16 and for the function 
%    ff = chebfun(@(zk) 1/g(zk), [3.4,3.5]);
%    we see that only 18 degree polynomials is needed and the roots agrees
%    with the Fredholm techniques.
% 3. For the number of holes = 24 and for the function 
%    ff = chebfun(@(zk) 1/g(zk), [3.4,3.5]);
%    we see that only 13 degree polynomials is needed and the roots agrees
%    with the Fredholm techniques, while for the Fredholm determinant we need 
%    15 degree polynomial. It appears that the optimization function may
%    be needing lesser degree interpolant as compared to the Fredholm as
%    the number of holes grows.