function [zer0, val0, varargout] = zeros_fred_opti(func, method, chnkr, opts)
% zero_fred_opti returns zeros and the value of the objective function at
% 'zer0'.
% 
% Input Arguments are
%   * 'func' for objective function
%       ** 'f', 'fred', 'fredholm_determinant'  : for Fredholm determinant
%       ** 'ir','inv_res','inverse_resolvent'   : for Inverse of the Resolvent (1/<u_i,inv(A)v_i>)
%
%   * 'method' for root finding algorithm 
%       ** 'c', 'cheb', 'chebyshev'             : for Chebyshev interpolant
%       ** 'cm', 'comp', 'complex_muller'       : for Complex Muller
%       ** 'n', 'newt', 'newton'                : for Newton
% 
%   * 'chnkr' geometric object for which you are computing the eigenvalues.
%
% Optional input arguments:
%   opts  - options structure. available options (default settings)
%       * [opts.cheb1, opts.cheb1] = Interval where we are interpolating 
%                                    the objective function using Chebfun
%       * opts.freddeg = Degree of Chebyshev Interpolant for the Fredholm determinant
%       * opts.irdeg = Degree of Chebyshev Interpolant for the Inverse of the Resolvent
%       * opts.cmguess = Initial Guess for Complex Muller
%       * opts.cmiter = Iterations for Complex Muller
%       * opts.newtguess = Initial Guess for Newton
%       * opts.newtiter = Iterations for Newton
%       * opts.u = u(=[u_i]) in the Inverse of the Resolvent (1/<u_i,inv(A)v_i>)
%       * opts.v = v(=[v_i]) in the Inverse of the Resolvent (1/<u_i,inv(A)v_i>)
%       * opts.nd = Number of u's and v's
%       * opts.speed = 'slow' for Slow evaluation and 'fast' for Fast evaluation






cheb1 = 3;                  % Left point of the interval where chebyshev expansion is done
cheb2 = 4;                  % Right point of the interval where chebyshev expansion is done
freddeg = 16;               % Degree of Chebyshev Interpolant for the Fredholm determinant
irdeg = 35;                 % Degree of Chebyshev Interpolant for the Inverse Resolvent
cmguess = [3 3.1 3.5];      % Initial Guess for Complex Muller
cmiter = 10;                % Iterations for Complex Muller
newtguess = 3;              % Initial Guess for Newton
newtiter = 10;              % Iterations for Newton
nd = 1;                     % Number of u's and v's
speed = 'fast';

if isfield(opts,'cheb1')
    cheb1 = opts.cheb1;
end
if isfield(opts,'cheb2')
    cheb2 = opts.cheb2;
end
if isfield(opts,'freddeg')
    freddeg = opts.freddeg;
end
if isfield(opts,'irdeg')
    irdeg = opts.irdeg;
end
if isfield(opts,'cmguess')
    cmguess = opts.cmguess;
end
if isfield(opts,'newtguess')
    newtguess = opts.newtguess;
end
if isfield(opts,'newtiter')
    newtiter = opts.newtiter;
end
if isfield(opts,'cmiter')
    cmiter = opts.cmiter;
end
if isfield(opts,'nd')
    nd = opts.nd;
end
xx = chnkr.r(1,:);
yy = chnkr.r(2,:);
if isfield(opts,'u')
    u = opts.u;
else
    u = zeros(nd, length(xx));
    v = zeros(nd, length(xx));
    for ii = 1:nd
        u(ii,:) = (xx.^(2*(1/(ii)) + 2) - yy.^(2*(1/(ii)) + 2) + yy.^(2*(1/(ii)) + 3)).*sin(xx) + cos(xx).*(0.5*yy.^(2*(1/(ii)) + 2) + 0.3*yy);
        v(ii,:) = (xx.^(2*(1/(ii)) + 3) - yy.^(3*(1/(ii)) + 2) + cos(xx.^(2*(1/(ii)) + 2))).*cos(-yy) + sin(yy).*(xx.^(2*(1/(ii)) + 2) + xx);
    end
end
if isfield(opts,'v')
    v = opts.v;
else
    [nd,~] = size(u);
    v = zeros(nd, length(xx));
    for ii = 1:nd
        v(ii,:) = (xx.^(2*(1/(ii)) + 3) - yy.^(3*(1/(ii)) + 2) + cos(xx.^(2*(1/(ii)) + 2))).*cos(-yy) + sin(yy).*(xx.^(2*(1/(ii)) + 2) + xx);
    end
end
if isfield(opts,'speed')
    speed = opts.speed;
end

obj = [];
obj.u = u;
obj.v = v;
obj.speed = speed;

switch lower(method)
    case {'c', 'cheb', 'chebyshev'}
        switch lower(func)
            case {'f', 'fred', 'fredholm_determinant'}
                f = chebfun(@(zk) fred_opti(zk, chnkr, func, obj), [cheb1, cheb2], freddeg);
            case{'ir','inverse_resolvent', 'inv_res'}
                f = chebfun(@(zk) fred_opti(zk, chnkr, func, obj), [cheb1, cheb2], irdeg);
        end
        % change the coding so that the root with lowest imaginary part is
        % only returned
        rts = roots(f, 'complex');
        imag_parts = abs(imag(rts));
        [~, idx] = sort(imag_parts, 'ascend');
        sorted_rts = rts(idx);
        zer0 = sorted_rts(1);
        val0 = f(zer0);
        if nargout == 3
           varargout{1} = f;
        elseif nargout == 4
           varargout{1} = f;
           varargout{2} = rts;
        end
    case {'cm', 'comp', 'complex_muller'}
        x0 = cmguess(1); y0 = fred_opti(x0, chnkr, func, obj); 
        x1 = cmguess(2); y1 = fred_opti(x1, chnkr, func, obj); 
        x2 = cmguess(3); y2 = fred_opti(x2, chnkr, func, obj); 
        cmz = zeros(1, cmiter);
        cmy = zeros(1, cmiter);
        if (abs(y0) < 1e-12)
            zer0 = x0;
            val0 = y0;
        elseif (abs(y1) < 1e-12)
            zer0 = x1;
            val0 = y1;
        elseif (abs(y2) < 1e-12)
            zer0 = x2;
            val0 = y2;
        else 
           for ii = 1:cmiter
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
               x2 = x3; y2 = fred_opti(x2, chnkr, func, obj); 
               cmz(ii) = x2; cmy(ii) = y2;
               if (abs(y2) < 1e-12)
                    break;
               end
           end
           zer0 = x2;
           val0 = y2;
           if nargout == 3    %%% MANAS: How can we ensure if people dont ask for third entry but ask for fourth, then we only print the things that have been asked for  
                varargout{1} = ii;
           end
           if nargout == 4
               varargout{1} = ii;
               varargout{2} = [cmz.' cmy.'];
           end
        end

    case {'n', 'newt', 'newton'}
        val = fred_opti(newtguess, chnkr, func, obj);
        if abs(val)<1e-12
            zer0 = newtguess;
            val0 = val;
        else
            newtz = zeros(1, newtiter);
            newty = zeros(1, newtiter);
            for ii = 1:newtiter
               [g, derg] =  fred_opti(newtguess, chnkr, func, obj);
               newtguess = newtguess - g/derg;
               val = g;
               newtz(ii) = newtguess;
               newty(ii) = val;
               if (abs(val) < 1e-12)
                    break;
              end
            end
            zer0 = newtguess;
            val0 = val;
            if nargout == 3 
                varargout{1} = ii;
            end
            if nargout == 4
                varargout{1} = ii;
                varargout{2} = [newtz.' newty.'];
            end
        end
    otherwise
         fprintf('Input is invalid');
end

end