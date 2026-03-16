%%%%%%%%%%  SLOW CODE FOR FREDHOLM DETERMINANT USING CHEBYSHEV %%%%%%%%%
clear;
chnkreps = get_chunker(2);
% ks0 = tic;
% g = @(zk) get_determinant(zk, chnkreps);  
% f1 = get_chebyshev(g, 3.3, 3.5, 16);
% %%% problem with get_chebyshev is that it's coefficients are coefficients
% %%% from Chebyshev basis and hence we are unable to find roots.
% tks = toc(ks0);
t1 = tic;
f = chebfun(@(zk) get_determinant(zk, chnkreps), [3.3,3.5], 16);
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%   CHEBYSHEV NODES  %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ns = 1:5:200;
% err = zeros(size(ns));
% for jj = 1:length(ns)    
%     n = ns(jj);
%     t = 0:1:n;
%     x = cos(2*pi*t/(2*n + 1));              % Interpolating points
%     for kk = 1:length(t)
%         g(kk) = f(x(kk));                   % g = f(cos(t))
%     end
%     A = cos( (2 * pi / (2 * n + 1)) * [0:1:n].' * [0:1:n] );
%     d = 1/(2*n+1) * A * (g.*[1 2*ones(1,n)]).';
%     p = @(x) (cos(x*[0:1:n]).*[1 2*ones(1,n)])*d;
%     k = 10;
%     w = 0:pi/k:pi;
%     for ii = 1:length(w)
%         z(ii) = p(w(ii));
%         y(ii) = f(cos(w(ii)));
%     end
%     err(jj) = abs(max(y - z));
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  FREDHOLM DETERMINANT  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% function det1 = get_determinant(zk, chnkr1)
% 
%     Dk = 2*kernel('helm', 'd', zk);  
%     A = chunkermat(chnkr1, Dk);   
% %%
%     A = A + eye(chnkr1.npt);
% %%
% 
%     det1 = det(A);
% end