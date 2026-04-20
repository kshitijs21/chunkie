% f = @(x) sin(x);
% g = @(x) cos(x);
% maybe m is not needed
% m = 2;
chnkr = get_chunker(2);
xx = chnkr.r(1,:);
yy = chnkr.r(2,:);
u = (xx.^2 - yy.^2 + yy).*sin(xx) + cos(xx).*(0.5*yy.^2 + 0.3*yy);
v = (xx - yy.^3 + xx.^2).*cos(-yy) + sin(yy).*(xx.^2 + xx);
x0 = 3.449274696865717; x1 = 3.4479274686865717; x2 = 3.449674686865787;

opts = [];
opts.cmguess = [x0 x1 x2];
opts.u = u;
opts.v = v;
opts.cmiter = 50;
t1 = tic;
[zk0, val, niter, zzz] = zeros_fred_opti('ir', 'cm', chnkr, opts);
telg = toc(t1);
fprintf('Time taken by One above all = %d\n', telg);


A = @(zk) get_matrix(zk, chnkr);
g = @(zk) ones(1, length(u))*( ( (u.').*( A(zk)\(v.') ) ).*(chnkr.wts(:)) );
ff = @(zk) 1/g(zk);

x0 = 3.409274696865717; x1 = 3.409274686865717; x2 = 3.409274686865787;
y0 = ff(x0); y1 = ff(x1); y2 = ff(x2);
n = 50; % number of iterations
zzo = zeros(1,n);
yyo = zeros(1,n);
t2 = tic;
for i = 1:n
    h0 = x1 - x0; 
    h1 = x2 - x1;
    delta0 = (y1 - y0)/h0;
    delta1 = (y2 - y1)/h1;
    a = (delta1 - delta0)/(h1 + h0); 
    b = a*h1 + delta1;
    c = y2;

    x3 = x2 - 2*c/(b + sign(b)*sqrt(b*b - 4*a*c));

    x0 = x1; y0 = y1;
    x1 = x2; y1 = y2;
    x2 = x3; y2 = ff(x2);
    zzo(i) = x2; yyo(i) = y2; 
    if (abs(y2) < 1e-12)
        break;
    end
end
told = toc(t2);
niters = i;
fprintf('Time taken by old coding = %d\n', told);

%%%%% older method is because of the slow code and one above all is due to
%%%%% the fast code hence there is a slight difference.
 

function [A] = get_matrix(zk, chnkr1)

    Dk = 2*kernel('helm', 'd', zk);  
    A = chunkermat(chnkr1, Dk);   
%%
    A = A + eye(chnkr1.npt);
end
