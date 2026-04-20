%%%%%%%%%%  SLOW CODE FOR FREDHOLM DETERMINANT USING CHEBYSHEV %%%%%%%%%
%clear;
ll = [2 2.6 3];
% ll = 3;
del = [0.05 0.1 0.2];
zer0 = 3.235377405473126e+00;   % root corresponding to l = 1.5
% zer0 = 3.396302810423935; % root corresponding to l = 2.6
niters = zeros(length(ll), length(del));
zero0 = zeros(length(ll), length(del));
for ii = 1:length(ll)
    [chnkr, nh] = get_chunker(ll(ii));
    f = @(zk) get_determinant(zk, chnkr);
    fprintf('Hole = %d:\n', nh);
 
    for jj = 1:length(del)
        kk = del(jj);
        z1 = zer0; bb = z1;
        z0 = z1-kk; aa = z0; 
        z2 = z1+kk; cc = z2;
        y0 = f(z0); y1 = f(z1); y2 = f(z2);  
        nit = 0;
        t1 = tic;
        if (abs(y0) < 1e-12)
            zer0 = z0;
        elseif (abs(y1) < 1e-12)
            zer0 = z1;
        elseif (abs(y2) < 1e-12)
            zer0 = z2;
        else 
            
            while(abs(y2) > 1e-9)
                nit = nit+1;
                h0 = z1 - z0; 
                h1 = z2 - z1;
                delta0 = (y1 - y0)/h0;
                delta1 = (y2 - y1)/h1;
                a = (delta1 - delta0)/(h1 + h0); 
                b = a*h1 + delta1;
                c = y2;

                z3 = z2 - 2*c/(b + sign(b)*sqrt(b*b - 4*a*c));

                z0 = z1; y0 = y1;
                z1 = z2; y1 = y2;
                z2 = z3; y2 = f(z2);   % 1 function evaluation for each i 
            end
            niters(ii,jj) = nit;
            zero0(ii,jj) = z2;
            y1 = y2;
        end
        timE = toc(t1);
        save(['Parameter_Complex_Muller_with_' num2str(nh) '_holes_and_initial_guess_(' num2str(aa) ', ' num2str(bb) ', ' num2str(cc) ').mat'], 'z2', 'y2','f', 'nit', 'timE');   
    
        fprintf('  del = %d, Root = %d, fn value = %d, iterations = %d, # fn eval = %d, time = %d \n', kk, z2, y2, nit, nit +3, timE);
    end
    zer0 = z2;

end

function det1 = get_determinant(zk, chnkr1)

    Dk = 2*kernel('helm', 'd', zk);  
    A = chunkermat(chnkr1, Dk);   

    A = A + eye(chnkr1.npt);

    det1 = det(A);
end


