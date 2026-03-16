
% Interpolating the optimization function requires around 300 degree
% polynomials
% ncheb = 30;
chnkr = get_chunker(2);
xx = chnkr.r(1,:);
yy = chnkr.r(2,:);
u = (xx.^2 - yy.^2 + yy).*sin(xx) + cos(xx).*(0.5*yy.^2 + 0.3*yy);
v = (xx - yy.^3 + xx.^2).*cos(-yy) + sin(yy).*(xx.^2 + xx);
A = @(zk) get_matrix(zk, chnkr);
g = @(zk) ones(1, length(u))*( ( (u.').*( A(zk)\(v.') ) ).*(chnkr.wts(:)) );
ff = chebfun(@(zk) 1/g(zk), [3,3.5]);
%%
plot(imag(ff), 'k.');
rts = roots(ff, 'complex');

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