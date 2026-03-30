%%%%%%%%%%%%%%%%     SHITTIEST CODE    %%%%%%%%%%%%%%%%%%%%%%%%%


function varargout = parameter_fred_opti(fn, z0, chnkr, method, opts)

% parameter_fred_opti returns parameters required corresponding to the 
% root finding algorithm of the objective function which are 
%   * Chebyshev
%       ** Interval of interpolant
%       ** Degree of interpolant
%   * Complex Muller
%       ** Three initial guesses
% 
% Input Arguments are
%   * 'fn' for objective function
%       ** 'f', 'fred', 'fredholm_determinant'  : for Fredholm determinant
%       ** 'ir','inv_res','inverse_resolvent'   : for Inverse of the Resolvent (1/<u_i,inv(A)v_i>)
%
%   * 'z0' for initial guess
% 
%   * 'chnkr' geometric object on which parameters is needed.
%
%   * 'method' for root finding algorithm 
%       ** 'c', 'cheb', 'chebyshev'             : for Chebyshev interpolant
%       ** 'cm', 'comp', 'complex_muller'       : for Complex Muller
%
% Optional input arguments:
%   opts  - options structure. available options (default settings)
%       * opts.tol_in   =  Tolerance in the function such that we have
%                          $max(f(zz) < tol_in$
%       * opts.tol_deg  =  Tolerance for the degree of Chebyshev Interpolant 
%       * opts.tol_cm   =  Tolerance for the initial guess via Complex Muller 
%       * opts.delta    =  Initial guess to how large the interval you want
%       * opts.deg      =  Guess of Degree of interpolant 
%       * opts.cmguess  =  Initial Guess for Complex Muller
%
% Output Arguments are
%   * when 'method' is 'chebyshev'
%       ** deg      : Degree of the Chebyshev interpolant to be choosen
%       ** z0       : center of the interval [z0-delta, z0+delta]
%       ** delta    : radius of the interval

%
%   * when 'method' is 'complex_muller'
%       ** [z0-delta, z0, z0+delta] : Three initial guesses

%%%% We are using the arguments similar to the Bisection method for root
%%%% finding.


delta = 1;   
tol_in = 1e-7;              
tol_deg = 1e-15;
tol_cm = 1e-7;
deg = 3;
cmguess = [z0 z0-delta z0+delta];     


if isfield(opts,'tol_in')
    tol_in = opts.tol_in;
end
if isfield(opts,'tol_cm')
    tol_cm = opts.tol_cm;
end
if isfield(opts,'tol_deg')
    tol_deg = opts.tol_deg;
end
if isfield(opts,'delta')
    delta = opts.delta;
end
if isfield(opts,'deg')
    deg = opts.deg;
end
if isfield(opts,'cmguess')
    cmguess = opts.cmguess;
end

f = @(zk) fred_opti(zk, chnkr, fn);

switch lower(method)
    case {'c', 'cheb', 'chebyshev'}

    a = z0-delta; b = z0+delta;
    
    % Finding the intervals [a,b] for interpolation
    
    zz = [z0,a,b];
    fa = abs(f(zz(1))); fb = abs(f(zz(2))); f0 = abs(f(zz(2)));
    yy = [f0, fa, fb];
    [~, idx] = sort(yy, 'ascend');
    z0 = zz(idx(1)); 
    while (yy(1)>tol_in) || (yy(2)>tol_in) || (yy(3)>tol_in)
        delta = delta/2;
        a = z0-delta; b = z0+delta;
        fa = abs(f(a)); fb = abs(f(b)); f0 = abs(f(z0));
        yy = [f0, fa, fb];
        zz = [z0,a,b];
        [~, idx] = sort(yy, 'ascend');
        z0 = zz(idx(1));
    end
    

    del = delta/2;
    a = z0-del; b = z0+del;
    
    % Finding the degree of interpolant
    
    ff = @(x) f( (b-a)*x/2 + (b+a)/2 );
   
    t = 0:1:deg;
    x = cos(2*pi*t/(2*deg + 1));  
    gg = zeros(1,deg);
    for jj = 1:deg
        gg(jj) = ff(x(jj));
    end
    d_deg = 1/(2*deg+1) * (gg.*[1 2*ones(1,deg)]) * cos((1/(2*deg+1))*2*pi*deg*t.');
   
    while (abs(d_deg)>tol_deg)
        hh = ceil(log(abs(d_deg)/tol_deg));
        deg = deg + hh;
        ggw = zeros(1, deg);
        t = 0:1:deg;
        x = cos(2*pi*t/(2*deg + 1)); 
        for ll = 1:length(x)
            ggw(ll) = ff(x(ll));
        end
        d_deg = 1/(2*deg+1) * (ggw.*[1 2*ones(1,deg)]) * cos(2*pi*deg/(2*deg+1)*t.');
    end
    
    if narargout == 1
        varargout{1} = deg;

    elseif narargout == 2
        varargout{1} = deg;
        varargout{2} = z0;

    else
        varargout{1} = deg;
        varargout{2} = z0; 
        varargout{3} = del; 
    end
    
    case {'cm', 'comp', 'complex_muller'}
        
        fa = abs(f(cmguess(1))); fb = abs(f(cmguess(2))); f0 = abs(f(cmguess(3)));
        yy = [f0, fa, fb];
        [~, idx] = sort(yy, 'ascend');      
        z0 = cmguess(idx(1)); 
        while (yy(1)>tol_cm) || (yy(2)>tol_cm) || (yy(3)>tol_cm)
            delta = delta/2;
            a = z0-delta; b = z0+delta;
            zz = [z0,a,b];
            fa = abs(f(a)); fb = abs(f(b)); f0 = abs(f(z0));
            yy = [f0, fa, fb];
            [~, idx] = sort(yy, 'ascend');
            z0 = zz(idx(1)); 
        end
        del = delta;
        varargout = [z0-del,z0,z0+del];

    otherwise
        fprintf('Invalid method');
end
end