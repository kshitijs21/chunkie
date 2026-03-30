%%%%%%%%     CENTRAL DIFFERENCE TEST FOR DERIVATIVE OF DETERMINANT   %%%%%%%%
% clear;
% chnkr = get_chunker(1);
% det = @(zk) get_determinant(zk, chnkr);
% derdet = @(zk) get_derdeterminant(zk, chnkr);
% 
% zk = 1:0.2:50;
% err = zeros(1,length(zk));
% err1 = zeros(1,length(zk));
% h0 = 1e-4;
% zz = 3.2;
% for i = 1:length(zk)
%      h = h0/zk(i);
%      Dp = det(zz + h);  
%      Dm = det(zz - h);
%      derp = derdet(zz);
%      err(i) = norm((Dp - Dm)/h - 2*derp);
% end
% loglog(zk, err, 'k.');
% hold on; loglog(zk, 150*(h0./zk).^2, 'b.');

%%%%%%%%     NEWTON FOR DERIVATIVE OF DETERMINANT   %%%%%%%%

chnkreps = get_chunker(3);   % 1 function eval for chunker
% z0 = 3.4093; 
z0 = 3.44976372033;
n = 10; % number of iterations
for i = 1:n
    [z, y] = update_iterate(z0, chnkreps);   % 1 function eval for Newton which has 5 inherent function eval in update_iterate
    fprintf('i=%d  f=%d,  z0=%d\n',i, abs(y), real(z0));
    z0 = z;
    if (abs(y) < 1e-12)
        break;
    end
end
niters = i;
fprintf('Func eval=%d\n', 5*niters);


function [z, y] = update_iterate(zk, chnkr1)
    tdk = tic;
    Dk = 2*kernel('helm', 'd', zk);  % 1 function eval in update_iterate
    t1 = toc(tdk);
    fprintf('time for kernel =%d\n', t1);
    ta = tic;
    A = chunkermat(chnkr1, Dk);   % 2 function eval in update_iterate
    t2 = toc(ta);
    fprintf('time for chunkermat =%d\n', t2);
    A = A + eye(chnkr1.npt);
    
    y = det(A);     % 3 function eval in update_iterate
    
    tderdk = tic;
    derDk = 2*kernel('helm', 'fd_d', zk);   % 4 function eval in update_iterate
    t3 = toc(tdk);
    fprintf('time for derkernel =%d\n', t3);
    tdera = tic;
    derA = chunkermat(chnkr1, derDk);       % 5 function eval in update_iterate
    t4 = toc(tdera);
    fprintf('time for derchunkermat =%d\n', t4);
    z = zk - 1/trace(derA*inv(A));
end

% Newton needs more number of steps as the number of holes increases.

% Inside of update_iterate we are doing 5 iterations.

% Total number of function evaluations that we are doing is 5*niters.
