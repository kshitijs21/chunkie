function [f, varargout] = fred_opti(zk, chnkr, func, opts)
% fred_opti returns either the Fredholm determinant and its derivative or
% the optimisation function and its derivative based on the value of 'choice'.
%
% Input Arguments are
%   * 'zk' point where you are computing the function value.
%
%   * 'chnkr' geometric object over which you are computing the function value.
%
%   * 'func' is a character such that it takes the following values
%       ** 'f', 'fred', 'fredholm_determinant': for Fredholm determinant and it's derivative
%       ** 'ir','inv_res','inverse_resolvent' : for Inverse of the Resolvent (1/<u_i,inv(A)v_i>) and its derivative
% 
%   * 'opts' is a struct with the entries 'u' and 'v' such that 
%       ** opts.u is a (nd, chnkr.npt) matrix
%       ** opts.v is a (nd, chnkr.npt) matrix
%       ** opts.speed: 'slow' for Slow evaluation and 'fast' for Fast evaluation
%       ** opts.nd = Number of u's and v's
%
% Output variables are
%   * 'f' function value.
%
%   * 'varargout{1}' derivative of the function value.
%
%   * 'varargout{2}' time taken for function evaluation
% 
%   * 'varargout{3}' time taken for derivative of function evaluation


nd  = 1;
xx = chnkr.r(1,:);
yy = chnkr.r(2,:);

if isfield(opts,'nd')
    nd = opts.nd;
end
rng(1);


% HW. Resolve the function u and v properly on the boundary
if isfield(opts,'u')
    u = opts.u;
else
    u = zeros(nd, length(xx));
    for ii = 1:nd
        rr = ii + rand(6,1);
        rr2 = rand(3,1);
        
        uu = @(x, y) (rr(1)*x + rr(2)*y + rr(3)*(x.*y) + rr(4)*x.^2 + rr(5)*y.^2 + rr(6)) .* sin(pi*(rr2(1)*x  +rr2(2)*y + rr2(3))) ;
        u(ii,:) = uu(xx,yy);
     end
end

if isfield(opts,'v')
    v = opts.v;
else
    [nd,~] = size(u);
    v = zeros(nd, length(xx));
    for ii = 1:nd
        rrr = nd + ii *rand(6,1);
        rr2 = rand(3,1);
        
        vv = @(x, y) (rrr(1)*x.^2 + rrr(2)*y.^3 + rrr(3)*(x.*y) + rrr(4)*x + rrr(5)*y.^2 + rrr(6)) .* sin(pi*(rr2(1)*x  +rr2(2)*y + rr2(3))) ;
        v(ii,:) = vv(xx,yy);
    end
end


speed = 'i_am_the_knightmayor';
if isfield(opts, 'speed')
    speed = opts.speed;
end
opts_flam = [];
opts_flam.flamtype = 'rskelf';
opts_flam.forceproxy = true;
opts_flam.occ = 200;

switch lower(func)
    case {'f', 'fred', 'fredholm_determinant'}
        Dk = 2*kernel('helm', 'd', zk);
        if strcmpi (speed, 'slow')  || (chnkr.npt < 4000 && ~strcmpi (speed, 'fast'))
            t1 = tic;
            A = chunkermat(chnkr, Dk);   
            A = A + eye(chnkr.npt);
            f = det(A);
            timef = toc(t1);
            if nargout >= 2
                t2 = tic;
                derDk = 2*kernel('helm', 'fd_d', zk);  
                derA = chunkermat(chnkr, derDk);    
                varargout{1} = f*trace(derA*inv(A));
                timeder = toc(t2);
                
                if nargout == 3
                    varargout{2} = timef;
               
                elseif nargout == 4
                    varargout{2} = timef;
                    varargout{3} = timeder;
                else
                    varargout{2} = timef;
                    varargout{3} = timeder;
                    fprintf('Only four output are possible\n');
                end
            end

        else
                
            dval = 1;
            t1 = tic;
            F = chunkerflam(chnkr, Dk, dval, opts_flam);
            f = exp(rskelf_logdet(F));
            timef = toc(t1);
            
            if nargout >= 2
                t2 = tic;
                A = chunkermat(chnkr, Dk);   
                A = A + eye(chnkr.npt);
                derDk = 2*kernel('helm', 'fd_d', zk);  
                derA = chunkermat(chnkr, derDk);    
                varargout{1} = f*trace(derA*inv(A));
                 timeder = toc(t2);
                 
                if nargout == 3
                    varargout{2} = timef;
               
                elseif nargout == 4
                    varargout{2} = timef;
                    varargout{3} = timeder;
                else
                    varargout{2} = timef;
                    varargout{3} = timeder;
                    fprintf('Only four output are possible\n');
                end
            end
        end
      
    case{'ir', 'inv_res', 'inverse_resolvent'}
        [ru, cu] = size(u);
        [rv, cv] = size(v);
        if cu~=chnkr.npt || cv~=chnkr.npt || cu~=cv || ru~=rv
            fprintf('Dimensions of u and v are not correct\n');
        else
            Dk = 2*kernel('helm', 'd', zk); 

            if strcmpi (speed, 'slow') || (chnkr.npt < 4000 && ~strcmpi (speed, 'fast'))
                
                t1 = tic;
                A = chunkermat(chnkr, Dk);    
                A = A + eye(chnkr.npt);
                opti = sum( ones(1, length(u))*( ( (u.').*( A\(v.') ) ).*(chnkr.wts(:)) ) );
                f = 1/opti;
                timef = toc(t1);
                
                if nargout >=2 
                    t2 = tic;
                    Derk = 2*kernel('helmholtz','fd_d',zk);
                    derA = chunkermat(chnkr, Derk);
                    deropti = ones(1, length(u))*...
                                (   (u.').*(  A \( derA*(A\(v.')) )  ).*(chnkr.wts(:))   );
                    varargout{1} = deropti/(opti^2);
                    timeder = toc(t2);
                    
                    if nargout == 3
                    varargout{2} = timef;
               
                    elseif nargout == 4
                    varargout{2} = timef;
                    varargout{3} = timeder;
                    elseif nargout > 4 
                    varargout{2} = timef;
                    varargout{3} = timeder;
                    fprintf('fred_opti returns four outputs\n');
                    end
                end
            

            else
                t1 = tic;
                F = chunkerflam(chnkr, Dk, 1.0, opts_flam); 
                sol = (rskelf_sv(F,v.'));
                opp = (u.').*sol.*chnkr.wts(:);
                opti = sum(opp(:));  
                f = 1/opti;
                timef = toc(t1);
                
                if nargout >= 2
                    t2 = tic;
                    Derk = 2*kernel('helmholtz','fd_d',zk);
                    choice0 = [];
                    choice0.corrections = true;
                    cormat = chunkermat(chnkr, Derk, choice0);
                    choice = [];
                    choice.forcesmooth = true;
                    choice.cormat = cormat;
                    u_eval = zeros(chnkr.npt, ru);
                    for ii = 1:ru
                        u_eval(:,ii) = chunkerkerneval(chnkr,Derk,sol(:,ii),chnkr,choice);
                    end
                    solder = rskelf_sv(F, u_eval); 
                    deropp = (u.').*(solder.*chnkr.wts(:));
                    deropti = sum(deropp(:));
                    varargout{1} = deropti/(opti^2);
                    timeder = toc(t2);
                    
                    if nargout == 3
                        varargout{2} = timef;
               
                    elseif nargout == 4
                        varargout{2} = timef;
                        varargout{3} = timeder;
                    elseif nargout > 4
                        varargout{2} = timef;
                        varargout{3} = timeder;
                        fprintf('fred_opti returns four outputs\n');
                    end
                end
            end
        end
    otherwise
        fprintf('Choice is invalid\n');
end
end