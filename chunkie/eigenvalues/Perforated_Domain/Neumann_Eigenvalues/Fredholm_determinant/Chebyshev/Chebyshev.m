%%%%%%%%%%  SLOW CODE FOR FREDHOLM DETERMINANT USING CHEBYSHEV %%%%%%%%%
clear;
chnkreps = get_chunker(2);

z0 = 3.4;
del = 0.1;
a = z0-del; b = z0 + del;
freddeg = 16;

t1 = tic;
f = chebfun(@(zk) get_determinant(zk, chnkreps), [a,b], freddeg);
timefuncdef = toc(t1);
t2 = tic;
plot(real(f), 'k.');
timefuncplot = toc(t2);
fprintf('Time in func def = %d, Time in func real plot = %d, ', timefuncdef, timefuncplot);
t3 = tic;
rts1 = roots(f, 'complex');
timerts = toc(t3);
fprintf('Time taken to find roots = %d', timerts);


function det1 = get_determinant(zk, chnkr1)

    Dk = 2*kernel('helm', 'd', zk);  
    A = chunkermat(chnkr1, Dk);   

    A = A + eye(chnkr1.npt);

    det1 = det(A);
end