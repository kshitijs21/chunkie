function [zk0, fval0, varargout] = zeros_fred_opti(funchoice, zerochoice, chnkr, obj)
% zero_fred_opti returns zeros and the value of the objective function at 'zero' d its derivative or
% the optimisation function and its derivative based. To get the Fredholm
% determinant
g = @(zk) fred_opti(zk, chnkr, funchoice, obj);
derg = @(zk) der_fred_opti(zk, chnkr, funchoice, obj);
U = [];

% g = fred_opti(zk, chnkr, funchoice, obj);
% [g, derg] = fred_opti(zk, chnkr, funchoice, obj);

if strcmp(zerochoice, 'c')
    if strcmp(funchoice, 'F')
        f = chebfun(@(zk) g(zk), [3.3,3.5], 16);
    else
        f = chebfun(@(zk) g(zk), [3.3,3.5], 50);
    end
    U.zero = roots(f, 'complex');
    U.val = f(U.zero);
    if nargout == 3
        varargout{1} = f;
    end
    
elseif strcmp(zerochoice, 'cm')
    x = input('Enter three initial guesses:\n');
    x0 = x(1); y0 = g(x(1)); 
    x1 = x(2); y1 = g(x(2)); 
    x2 = x(3); y2 = g(x(3));     % 3 function evaluation for each j
    if (abs(y0) < 1e-12)
        U.zero = x0;
        U.val = y0;
    elseif (abs(y1) < 1e-12)
        U.zero = x1;
        U.val = y1;
    elseif (abs(y2) < 1e-12)
        U.zero = x2;
        U.val = y2;
    else 
        n = input('Enter the number of iterations you want to do:\n');
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
            x2 = x3; y2 = g(x2);   % 1 function evaluation for each i 
            if (abs(y2) < 1e-12)
                U.zero = x2;
                U.val = y2;
                varargout{1} = i;
                break;
            end
        end
    end
elseif strcmp(zerochoice, 'n')
    z0 = input('Enter the initial guess:\n');
    val = g(z0);
    if abs(val)<1e-12
        U.zero = z0;
        U.val = val;
    else
        n = input('Enter the number of iterations you want:\n');
        for i = 1:n
            z0 = z0 - g(z0)/derg(z0);
            val = g(z0);
            if (abs(val) < 1e-12)
                U.zero = z0;
                U.val = val;
                varargout{1} = i;
                break;
            end
        end
    end
else
    fprintf('Input is invalid');
end
end