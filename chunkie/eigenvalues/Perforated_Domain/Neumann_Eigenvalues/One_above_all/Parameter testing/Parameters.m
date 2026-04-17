clear;
func = 'fred';
method = 'cm';
opts = [];


ll = [2 2.6 3];
% nh = [12 16 24];
del = [0.05 0.1 0.2];
% z00 = 3.44;         % replace by the roots of the previous value of ll, i.e., z0 = root of ll(i-1)
z00 = 3.235377405;
switch lower(method)
    case {'c', 'cheb', 'chebyshev'}
        for ii = 1:length(ll)
            [chnkr, nh] = get_chunker(ll(ii));
            fprintf('Hole = %d: \n', nh);
            for jj = 1:length(del)
              kk = del(jj);
              a = z00-kk;b = z00+kk;
                    
             f = chebfun(@(zk) fred_opti(zk, chnkr, func, opts), [a,b]);
                    
             rts = roots(f, 'complex');
             [nr,cr] = size(rts);
             aaa = 0;
             if (nr>=1) && (cr>=1)
                    aaa = 1;
                    imag_parts = abs(imag(rts));
                    [~, idx] = sort(imag_parts, 'ascend');
                    sorted_rts = rts(idx);
                   zer0 = sorted_rts(1,1);
                    
                   save(['Parameter_Cheb_' char(func) '_with_' num2str(nh) '_holes_on_interval_[' num2str(a) ',' num2str(b) '].mat'], 'rts', 'zer0','f');
            
             end
             if aaa == 0
                    fprintf('   No roots in the interval [%d,%d]\n', a, b);
                    save(['No_roots_Cheb_' char(func) '_with_' num2str(nh) '_holes_on_interval_[' num2str(a) ',' num2str(b) '].mat'], 'f');
                    zer0 = z00;
             end
            end
            z00 = zer0;
        end
    case{'cm', 'comp', 'complex_muller'}
        niters = zeros(length(ll), length(del));       % Collection of number of iterations for holes in {12,16,24} and delta in {0.05,0.1,0.2}
        zero0 = zeros(length(ll), length(del));        % Collection of roots for holes in {12,16,24} and delta in {0.05,0.1,0.2}
        fun0 = zeros(length(ll), length(del));         % Collection of function value at roots for holes in {12,16,24} and delta in {0.05,0.1,0.2}
        nholes = zeros(length(ll), length(del));       % Collection of holes for holes in {12,16,24} and delta in {0.05,0.1,0.2}
        timE = zeros(length(ll), length(del));         % Collection of time taken by each iterations for holes in {12,16,24} and delta in {0.05,0.1,0.2}
        zit = zeros(length(ll), length(del), 50);      % Collection of iterates for holes in {12,16,24} and delta in {0.05,0.1,0.2}
        yit = zeros(length(ll), length(del), 50);      % Collection of function value at iterates for holes in {12,16,24} and delta in {0.05,0.1,0.2}

        for ii = 1:length(ll)
            [chnkr, nh] = get_chunker(ll(ii));
            f = @(zk) fred_opti(zk, chnkr, func, opts);
            fprintf('Hole = %d:\n', nh);
 
            for jj = 1:length(del)
               kk = del(jj);
               z1 = z00; bb = z1;
               z0 = z1-kk; aa = z0; 
               z2 = z1+kk; cc = z2; 
               y0 = f(z0); y1 = f(z1); y2 = f(z2);  
               nit = 0;
               t1 = tic;
               if (abs(y0) < 1e-12)
                   nit = 1;
                   z2 = z0;
               elseif (abs(y1) < 1e-12)
                   nit = 1; 
                   z2 = z1;
               elseif (abs(y2) < 1e-12)
                   nit = 1; 
                   z2 = z2;
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
                     zit(ii,jj,nit) = z2;
                     yit(ii,jj,nit) = y2;
                    end
                niters(ii,jj) = nit;
                zero0(ii,jj) = z2;
                fun0(ii,jj) = y2;
                zzit = zit(ii,jj,1:1:nit);
                yyit = yit(ii,jj,1:1:nit);

                y1 = y2;
               end
               tt = toc(t1);
               timE(ii,jj) = tt;
               nholes(ii,jj) = nh;
               save(['Parameter_Complex_Muller_' char(func) '_iterates_with_' num2str(nh) '_holes_and_initial_guess_(' num2str(aa) ', ' num2str(bb) ', ' num2str(cc) ').mat'], 'z2', 'y2','f', 'nit', 'timE', 'zzit', "yyit");   
               fprintf('  del = %d, Root = %d, fn value = %d, iterations = %d, # fn eval = %d, time = %d \n', kk, z2, y2, nit, nit +3, tt);
            end
            z00 = z2;

        end
        save(['Parameter_Complex_Muller_' char(func) '_ultra_saver.mat'], 'zero0', 'fun0', 'niters', 'timE', 'nholes', 'zit', 'yit');   
    otherwise
        fprintf('Invalid inputs\n');
end


%%%%%%% To get iterates in Complex Muller with hole corresponding to 'ii' hole and
%%%%%%% delta corresponding to 'del(jj)', write the following 
%%%%%%% zw = zit(ii,jj,:); 
%%%%%%% zz = zw(1,:);
%%%%%%% 'zz' hence gives the iterates corresponding to "l(ii)" and delta
%%%%%%% "del(jj)" and for the function values corresponding to
%%%%%%% these iterates just replace 'zit' with 'yit'.