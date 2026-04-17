clear;
func = 'ir';
method = 'c';
opts = [];
opts.nd = 1;

fprintf('Objective funtion = %s\n', func);

ll = [2 2.6 3];
% nh = [12 16 24];
del = [0.05 0.1 0.2];
z00 = 3.235377405;
% z00 = 3.44976;
% z00 = 3.441582096927;
% z00 = 3.4497637315;
switch lower(method)
    case {'c', 'cheb', 'chebyshev'}
        fprintf('Root Finding Algorithm = Chebyshev \n');
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
              save([ 'Parameter_new_Chebysheff_' char(func) '_with_' num2str(nh) '_holes_on_interval_[' num2str(a) ',' num2str(b) '].mat'], 'rts', 'zero0','val0','f');
              
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
        fprintf('Root Finding Algorithm = Complex Muller\n');
        nholes = zeros(length(ll), length(del));       % Collection of holes for holes in {12,16,24} and delta in {0.05,0.1,0.2}
        timE = zeros(length(ll), length(del));         % Collection of time taken by each iterations for holes in {12,16,24} and delta in {0.05,0.1,0.2}
      
        for ii = 1:length(ll)
            [chnkr, nh] = get_chunker(ll(ii));
            fprintf('Hole = %d: Initial guess = %d \n', nh, z00);
 
            for jj = 1:length(del)
               delta = del(jj)
               z1 = z00; bb = z1;
               z0 = z1-delta; aa = z0; 
               z2 = z1-delta/2; cc = z2; 
               cguess = [z0 z2 z1];          % Assign initial guess here
               opts.cmguess = cguess;
               t1 = tic;
               [zero0, val0, nit, valit] = zeros_fred_opti(func, method, chnkr, opts);
               tt = toc(t1);
               timE(ii,jj) = tt;
               nholes(ii,jj) = nh;
               zzit = valit(1,:);
               yyit = valit(2,:);
               save(['Parameter_Complex_Muller_' char(func) '_iterates_with_' num2str(nh) '_holes_and_initial_guess_(' num2str(cguess(1)) ', ' num2str(cguess(2)) ', ' num2str(cguess(3)) ').mat'], 'zero0', 'val0', 'nit', 'timE', 'zzit', "yyit");   
               fprintf('  del = %d, Root = %d, fn value = %d, iterations = %d, #fn eval = %d, time = %d \n', delta, zero0, val0, nit, nit +3, tt);
            end
            z00 = zero0;

        end  

        case {'n', 'newt', 'newton'}
        fprintf('Root Finding Algorithm = Newton\n');
        % z00 = 3.44;
        nholes = zeros(1, length(ll));       % Collection of holes for holes in {12,16,24} and delta in {0.05,0.1,0.2}
        tim = zeros(1, length(ll));          % Collection of time taken by each iterations for holes in {12,16,24} and delta in {0.05,0.1,0.2}
        opts.newtiter = 50;
        
        for ii = 1:length(ll)
            [chnkr, nh] = get_chunker(ll(ii));
            fprintf('Hole = %d: Initial guess = %d \n', nh, z00);
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
               save(['Parameter_new_Newton_' char(func) '_iterates_with_' num2str(nh) '_holes_and_initial_guess_' num2str(aa) '.mat'], 'zero0', 'val0', 'nit', 'timE', 'zzit', "yyit", "derit");   
               fprintf('  Root = %d, fn value = %d, iterations = %d, time = %d \n', zero0, val0, nit, timE);
           
            z00 = zero0;

        end  

    case{'s', 'sec', 'secant'}
        fprintf('Root Finding Algorithm = Secant \n');
        nholes = zeros(length(ll), length(del));       % Collection of holes for holes in {12,16,24} and delta in {0.05,0.1,0.2}
        % timE = zeros(length(ll), length(del));         % Collection of time taken by each iterations for holes in {12,16,24} and delta in {0.05,0.1,0.2}
      
        for ii = 1:length(ll)
            [chnkr, nh] = get_chunker(ll(ii));
            fprintf('Hole = %d: Initial guess = %d \n', nh, z00);
 
            for jj = 1:length(del)
               delta = del(jj)
               z1 = z00; bb = z1;
               z0 = z1-delta; aa = z0; 
               z2 = z1+delta; cc = z2; 
               sguess = [z0 z1];             % Assign initial guess here
               opts.secguess = sguess;
               t1 = tic;
               [zero0, val0, nit, valit] = zeros_fred_opti(func, method, chnkr, opts);
               timE = toc(t1);
              
               nholes(ii,jj) = nh;
               zzit = valit(1,:);
               yyit = valit(2,:);
               save(['Parameter_new_Secant_' char(func) '_iterates_with_' num2str(nh) '_holes_and_initial_guess_(' num2str(sguess(1)) ', ' num2str(sguess(2)) ').mat'], 'zero0', 'val0', 'nit', 'timE', 'zzit', "yyit");   
               fprintf('  del = %d, Root = %d, fn value = %d, iterations = %d, # fn eval = %d, time = %d \n', delta, zero0, val0, nit, nit +2, timE);
            end
            z00 = zero0;

        end  

    otherwise
        fprintf('Invalid inputs\n');
end