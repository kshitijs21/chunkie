%%%%%%%%%%  SLOW CODE FOR FREDHOLM DETERMINANT USING CHEBYSHEV %%%%%%%%%
%clear;
ll = [2];
del = [0.05 0.1 0.2];
z0 = 3.235377405473126e+00;         % replace by the roots of the previous value of ll, i.e., z0 = root of ll(i-1)
for ii = 1:length(ll)
    [chnkr, nh] = get_chunker(ll(ii));
    fprintf('Hole = %d\n', nh);
    for jj = 1:length(del)
        kk = del(jj);
        a = z0-kk;b = z0+kk;
                    
        f = chebfun(@(zk) get_determinant(zk, chnkr), [a,b]);
                    
        rts = roots(f, 'complex');
        [nr,cr] = size(rts);
        aaa = 0;
        if (nr>=1) && (cr>=1)
            aaa = 1;
            imag_parts = abs(imag(rts));
            [~, idx] = sort(imag_parts, 'ascend');
            sorted_rts = rts(idx);
            zer0 = sorted_rts(1,1)
                    
            save(['Hyperparameter_Chebyshev_with_' num2str(nh) '_holes_on_interval_[' num2str(a) ',' num2str(b) '].mat'], 'rts', 'zer0','f');
            
        end
        if aaa == 0
            fprintf('   No roots in the interval [%d,%d]\n', a, b);
            save(['No_roots_Chebyshev_with_' num2str(nh) '_holes_on_interval_[' num2str(a) ',' num2str(b) '].mat'], 'f');
        end
    end
    z0 = zer0;
end

function det1 = get_determinant(zk, chnkr1)

    Dk = 2*kernel('helm', 'd', zk);  
    A = chunkermat(chnkr1, Dk);   

    A = A + eye(chnkr1.npt);

    det1 = det(A);
end


