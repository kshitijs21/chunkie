clear;
func = 'fred';
method = 's';
opts = [];



ll = [2.6 3];
% nh = [12 16 24];
del = [0.05 0.1 0.2];
% z00 = 3.235377405;
z00 = 3.44976;
% z00 = 3.441582096927;
switch lower(method)
    case {'c', 'cheb', 'chebyshev'}
        fprintf('Chebyshev \n');
        for ii = 1:length(ll)
            [chnkr, nh] = get_chunker(ll(ii));
            aa = 0; bb = 0;
            fprintf('Hole = %d: Initial guess = %d \n', nh, z00);
            for jj = 1:length(del)
              delta = del(jj)
              a = z00-delta;
              b = z00+delta;
              opts.cheb1 = a; opts.cheb2 = b;
              
              [zero0, val0, f, rts] = zeros_fred_opti(func, method, chnkr, opts);
              save([ char(opts.speed) '_Parameter_Chebysheff_' char(func) '_with_' num2str(nh) '_holes_on_interval_[' num2str(a) ',' num2str(b) '].mat'], 'rts', 'zero0','val0','f');
              
              [aa,bb] = size(rts);
              if (aa >=1) && (bb>=1)
                fprintf('  Interval = [%d, %d], delta = %d, zero = %d, deg = %d\n', a, b, delta, zero0, length(f)-1);   
                 z00 = zero0;
              else
                fprintf('  Interval = [%d, %d], delta = %d, zero = %s\n', a, b, delta, zero0);
              end
            
            end
         end
    case{'cm', 'comp', 'complex_muller'}
        fprintf('Complex Muller\n');
        nholes = zeros(length(ll), length(del));       % Collection of holes for holes in {12,16,24} and delta in {0.05,0.1,0.2}
        timE = zeros(length(ll), length(del));         % Collection of time taken by each iterations for holes in {12,16,24} and delta in {0.05,0.1,0.2}
      
        for ii = 1:length(ll)
            [chnkr, nh] = get_chunker(ll(ii));
            fprintf('Hole = %d:\n', nh);
 
            for jj = 1:length(del)
               delta = del(jj)
               z1 = z00; bb = z1;
               z0 = z1-delta; aa = z0; 
               z2 = z1+delta; cc = z2; 
               opts.cmguess = [z0 z1 z2];
               t1 = tic;
               [zero0, val0, nit, valit] = zeros_fred_opti(func, method, chnkr, opts);
               tt = toc(t1);
               timE(ii,jj) = tt;
               nholes(ii,jj) = nh;
               zzit = valit(1,:);
               yyit = valit(2,:);
               save([ char(opts.speed) '_Parameter_Complex_Muller_' char(func) '_iterates_with_' num2str(nh) '_holes_and_initial_guess_(' num2str(aa) ', ' num2str(bb) ', ' num2str(cc) ').mat'], 'zero0', 'val0', 'nit', 'timE', 'zzit', "yyit");   
               fprintf('  del = %d, Root = %d, fn value = %d, iterations = %d, # fn eval = %d, time = %d \n', delta, zero0, val0, nit, nit +3, tt);
            end
            z00 = zero0;

        end  

        case {'n', 'newt', 'newton'}
        fprintf('Newton method\n');
        z00 = 3.44;
        nholes = zeros(1, length(ll));       % Collection of holes for holes in {12,16,24} and delta in {0.05,0.1,0.2}
        tim = zeros(1, length(ll));         % Collection of time taken by each iterations for holes in {12,16,24} and delta in {0.05,0.1,0.2}
        opts.newtiter = 50;

        for ii = 1:length(ll)
            [chnkr, nh] = get_chunker(ll(ii));
            fprintf('Hole = %d:\n', nh);
            aa = z00;   
            opts.newtguess = aa;
               
               t1 = tic;
               [zero0, val0, nit, valit] = zeros_fred_opti(func, method, chnkr, opts);
               timE = toc(t1);
               tim(ii) = timE;
               nholes(ii) = nh;
               zzit = valit(1,:);
               yyit = valit(2,:);
               derit = valit(3,:);
               save(['Parameter_Newton_' char(func) '_iterates_with_' num2str(nh) '_holes_and_initial_guess_' num2str(aa) '.mat'], 'zero0', 'val0', 'nit', 'timE', 'zzit', "yyit", "derit");   
               fprintf('  Root = %d, fn value = %d, iterations = %d, time = %d \n', zero0, val0, nit, timE);
           
            z00 = zero0;

        end  

    case{'s', 'sec', 'secant'}
        fprintf('Secant method\n');
        nholes = zeros(length(ll), length(del));       % Collection of holes for holes in {12,16,24} and delta in {0.05,0.1,0.2}
        timE = zeros(length(ll), length(del));         % Collection of time taken by each iterations for holes in {12,16,24} and delta in {0.05,0.1,0.2}
      
        for ii = 1:length(ll)
            [chnkr, nh] = get_chunker(ll(ii));
            fprintf('Hole = %d:\n', nh);
 
            for jj = 1:length(del)
               delta = del(jj)
               z1 = z00; bb = z1;
               z0 = z1-delta; aa = z0; 
               z2 = z1+delta; cc = z2; 
               opts.secguess = [z0 z1];
               t1 = tic;
               [zero0, val0, nit, valit] = zeros_fred_opti(func, method, chnkr, opts);
               timE = toc(t1);
               tt(ii,jj) = timE;
               nholes(ii,jj) = nh;
               zzit = valit(1,:);
               yyit = valit(2,:);
               save(['Parameter_Secant_' char(func) '_iterates_with_' num2str(nh) '_holes_and_initial_guess_(' num2str(aa) ', ' num2str(cc) ').mat'], 'zero0', 'val0', 'nit', 'timE', 'zzit', "yyit");   
               fprintf('  del = %d, Root = %d, fn value = %d, iterations = %d, # fn eval = %d, time = %d \n', delta, zero0, val0, nit, nit +2, timE);
            end
            z00 = zero0;

        end  



    otherwise
        fprintf('Invalid inputs\n');
end