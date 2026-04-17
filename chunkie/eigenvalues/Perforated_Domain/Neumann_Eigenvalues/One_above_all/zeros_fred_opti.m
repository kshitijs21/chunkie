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
%       ** 's', 'sec', 'secant'                 : for Secant
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
%       * opts.secguess = Initial Guess for Secant
%       * opts.seciter = Iterations for Secant
%       * opts.u = u(=[u_i]) in the Inverse of the Resolvent (1/<u_i,inv(A)v_i>)
%       * opts.v = v(=[v_i]) in the Inverse of the Resolvent (1/<u_i,inv(A)v_i>)
%       * opts.nd = Number of u's and v's
%       * opts.speed = 'slow' for Slow evaluation and 'fast' for Fast evaluation






cheb1 = 3;                  % Left point of the interval where chebyshev expansion is done
cheb2 = 4;                  % Right point of the interval where chebyshev expansion is done
freddeg = 0;                % Degree of Chebyshev Interpolant for the Fredholm determinant
irdeg = 0;                  % Degree of Chebyshev Interpolant for the Inverse Resolvent
cmguess = [3 3.1 3.5];      % Initial Guess for Complex Muller
cmiter = 40;                % Iterations for Complex Muller
newtguess = 3;              % Initial Guess for Newton
newtiter = 10;              % Iterations for Newton
secguess = [3 3.1];         % Initial Guess for Secant method
seciter = 30;               % Iteration for Secant method
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
if isfield(opts,'secguess')
    secguess = opts.secguess;
end
if isfield(opts,'seciter')
    seciter = opts.seciter;
end
if isfield(opts,'cmiter')
    cmiter = opts.cmiter;
end
if isfield(opts,'nd')
    nd = opts.nd;
end
% xx = chnkr.r(1,:);
% yy = chnkr.r(2,:);
% if isfield(opts,'u')
%     u = opts.u;
% else
%     u = zeros(nd, length(xx));
%     v = zeros(nd, length(xx));
%     for ii = 1:nd
%         u(ii,:) = (xx.^(2*(1/(ii)) + 2) - yy.^(2*(1/(ii)) + 2) + yy.^(2*(1/(ii)) + 3)).*sin(xx) + cos(xx).*(0.5*yy.^(2*(1/(ii)) + 2) + 0.3*yy);
%         v(ii,:) = (xx.^(2*(1/(ii)) + 3) - yy.^(3*(1/(ii)) + 2) + cos(xx.^(2*(1/(ii)) + 2))).*cos(-yy) + sin(yy).*(xx.^(2*(1/(ii)) + 2) + xx);
%     end
% end

% if isfield(opts,'u')
%     u = opts.u;
% else
%     u = zeros(nd, length(xx));
%     for ii = 1:nd
%         rr = rand(length(xx),1);
% 
%         uu = @(x, y) (rr(1)*x + rr(2)*y + rr(3)*(x.*y) + rr(4)*x.^2 + rr(5)*y.^2 + rr(6)) .* sin(2*pi*x*rr(1)) ;
%         u(ii,:) = uu(xx,yy);
%         % vv = @(x, y) ((rrr.')*rrr + y*rr + (x.*(y.^3))*rrr + (x.^2)*rr + (x.^(2*ii))*rrr) * cos(y*rr + 2*pi*x*rr) ;
%         % v(ii,:) = vv(xx,yy);
%         % 
%         % u(ii,:) = (xx.^(2*(1/(ii)) + 2) - yy.^(2*(1/(ii)) + 2) + yy.^(2*(1/(ii)) + 3)).*sin(xx) + cos(xx).*(0.5*yy.^(2*(1/(ii)) + 2) + 0.3*yy);
%         % v(ii,:) = (xx.^(2*(1/(ii)) + 3) - yy.^(3*(1/(ii)) + 2) + cos(xx.^(2*(1/(ii)) + 2))).*cos(-yy) + sin(yy).*(xx.^(2*(1/(ii)) + 2) + xx);
%     end
% end
% 
% 
% if isfield(opts,'v')
%     v = opts.v;
% else
%     [nd,~] = size(u);
%     v = zeros(nd, length(xx));
%     for ii = 1:nd
% 
%         rrr = rand(length(xx),1);
%         vv = @(x, y) (rrr(1)*x.^2 + rrr(2)*y.^3 + rrr(3)*(x.*y) + rrr(4)*x + rrr(5)*y.^2 + rrr(6)) .* cos(2*pi*x*rrr(6)) ;
%         v(ii,:) = vv(xx,yy);
%     end
% end

if isfield(opts,'speed')
    speed = opts.speed;
end

obj = [];
if isfield(opts,'u')
obj.u = opts.u;
end
if isfield(opts,'v')
obj.v = opts.v;
end
obj.speed = speed;

switch lower(method)
    case {'c', 'cheb', 'chebyshev'}
        switch lower(func)
            case {'f', 'fred', 'fredholm_determinant'}
                if freddeg == 0
                    f = chebfun(@(zk) fred_opti(zk, chnkr, func, obj), [cheb1, cheb2]);
                else
                    f = chebfun(@(zk) fred_opti(zk, chnkr, func, obj), [cheb1, cheb2], freddeg);
                end
            case{'ir','inverse_resolvent', 'inv_res'}
                if irdeg == 0
                    f = chebfun(@(zk) fred_opti(zk, chnkr, func, obj), [cheb1, cheb2]);
                else
                    f = chebfun(@(zk) fred_opti(zk, chnkr, func, obj), [cheb1, cheb2], irdeg);
                end
        end
        % change the coding so that the root with lowest imaginary part is
        % only returned
        rts = roots(f, 'complex');
        [nr,cr] = size(rts);
        aaa = 0;
        if (nr>=1) && (cr>=1)
            aaa = 1;
            imag_parts = abs(imag(rts));
            [~, idx] = sort(imag_parts, 'ascend');
            sorted_rts = rts(idx);
            zer0 = sorted_rts(1);
            val0 = f(zer0);
        end
        if aaa == 0
            zer0 = 'No roots';
            val0 = 'No roots';
        end
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
            if nargout == 3    %%% MANAS: How can we ensure if people dont ask for third entry but ask for fourth, then we only print the things that have been asked for  
                varargout{1} = 0;
           end
           if nargout == 4
               varargout{1} = 0;
               varargout{2} = [cmz; cmy];
           end
        elseif (abs(y1) < 1e-12)
            zer0 = x1;
            val0 = y1;
            if nargout == 3    %%% MANAS: How can we ensure if people dont ask for third entry but ask for fourth, then we only print the things that have been asked for  
                varargout{1} = 0;
           end
           if nargout == 4
               varargout{1} = 0;
               varargout{2} = [cmz; cmy];
           end
        elseif (abs(y2) < 1e-12)
            zer0 = x2;
            val0 = y2;
            if nargout == 3    %%% MANAS: How can we ensure if people dont ask for third entry but ask for fourth, then we only print the things that have been asked for  
                varargout{1} = 0;
           end
           if nargout == 4
               varargout{1} = 0;
               varargout{2} = [cmz; cmy];
           end
        else 
            nit = 0;
            while(abs(y2) > 1e-9) && (cmiter >=nit)
               nit = nit+1;
           % for ii = 1:cmiter
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
               cmz(nit) = x2; cmy(nit) = y2;
               % if (abs(y2) < 1e-12)
               %      break;
               % end
           end
           zer0 = x2;
           val0 = y2;
           if nargout == 3    %%% MANAS: How can we ensure if people dont ask for third entry but ask for fourth, then we only print the things that have been asked for  
                varargout{1} = nit;
           end
           if nargout == 4
               varargout{1} = nit;
               varargout{2} = [cmz(1,1:1:nit); cmy(1,1:1:nit)];
           end
        end

    case {'n', 'newt', 'newton'}
        val = fred_opti(newtguess, chnkr, func, obj);
        newtz = zeros(1, newtiter);
        newtder = zeros(1, newtiter);
        newty = zeros(1, newtiter);
        if abs(val)<1e-12
            zer0 = newtguess;
            val0 = val;
            if nargout == 3 
                varargout{1} = 0;
            end
            newtz(1,1) = zer0;
            newty(1,1) = val0;
            if nargout == 4
                varargout{1} = 0;
                varargout{2} = [newtz;newty;newtder];
            end
        else
            nit = 0;
            while (abs(val) > 1e-9) && (newtiter >=nit)
                nit = nit+1;
                [g, derg] =  fred_opti(newtguess, chnkr, func, obj);
                newtguess = newtguess - g/derg;
                val = g;
                newtz(nit) = newtguess;
                newty(nit) = val;
                newtder(nit) = derg;
            end
            zer0 = newtguess;
            val0 = val;
            if nargout == 3 
                varargout{1} = nit;
            end
            if nargout == 4
                varargout{1} = nit;
                varargout{2} = [newtz; newty; newtder];
            end
        end
        %%%%% ADD secant method of root finding to this code.
    case {'s', 'sec', 'secant'}
        x0 = secguess(1); y0 = fred_opti(x0, chnkr, func, obj); 
        x1 = secguess(2); y1 = fred_opti(x1, chnkr, func, obj); 
        secz = zeros(1, seciter);
        secy = zeros(1, seciter);
        if (abs(y0) < 1e-12)
            zer0 = x0;
            val0 = y0;
            if nargout == 3    %%% MANAS: How can we ensure if people dont ask for third entry but ask for fourth, then we only print the things that have been asked for  
                varargout{1} = 0;
           end
           if nargout == 4
               varargout{1} = 0;
               varargout{2} = [secz; secy];
           end
        elseif (abs(y1) < 1e-12)
            zer0 = x1;
            val0 = y1;
            if nargout == 3    %%% MANAS: How can we ensure if people dont ask for third entry but ask for fourth, then we only print the things that have been asked for  
                varargout{1} = 0;
           end
           if nargout == 4
               varargout{1} = 0;
               varargout{2} = [secz; secy];
           end
        else
            nit = 0;
            while(abs(y1) > 1e-9) && (seciter >=nit)
               nit = nit+1;
               zer0 = x1 - y1*(x1-x0)/(y1-y0);
               val0 = fred_opti(zer0, chnkr, func, obj);
               x0 = x1;  y0 = y1;
               x1 = zer0;  y1 = val0;
               secz(nit) = x1; secy(nit) = y1;
            end
            if nargout == 3 
                varargout{1} = nit;
            end
            if nargout == 4
                varargout{1} = nit;
                varargout{2} = [secz; secy];
            end
        end
    otherwise
         fprintf('Input is invalid');
end

end