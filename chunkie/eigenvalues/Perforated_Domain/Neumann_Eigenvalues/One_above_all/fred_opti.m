function [f, varargout] = fred_opti(zk, chnkr, func, opts)
% fred_opti returns either the Fredholm determinant and its derivative or
% the optimisation function and its derivative based on the value of 'choice'.
%
% 'func' is a character such that it takes the following values
%       * 'f', 'fred', 'fredholm_determinant': for Fredholm determinant and it's derivative
%       * 'ir','inv_res','inverse_resolvent' : for Inverse of the Resolvent (1/<u_i,inv(A)v_i>) and its derivative
% 
% 'opts' is a struct with the entries 'u' and 'v' such that 
%       * opts.u is a (nd, chnkr.npt) matrix
%       * opts.v is a (nd, chnkr.npt) matrix
%       * opts.speed: 'slow' for Slow evaluation and 'fast' for Fast evaluation
%       * opts.nd = Number of u's and v's


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
        
        uu = @(x, y) (rr(1)*x + rr(2)*y + rr(3)*(x.*y) + rr(4)*x.^2 + rr(5)*y.^2 + rr(6)) .* sin(2*pi*x*rr(1)) ;
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
        
        vv = @(x, y) (rrr(1)*x.^2 + rrr(2)*y.^3 + rrr(3)*(x.*y) + rrr(4)*x + rrr(5)*y.^2 + rrr(6)) .* cos(2*pi*x*rrr(6)) ;
        v(ii,:) = vv(xx,yy);

        % rr = rand(length(xx),1);
        % rrr = rand(length(xx),1);
        % vv = @(x, y) ((rrr.')*rrr + y*rr + (x.*(y.^3))*rrr + (x.^2)*rr + (x.^(2*ii))*rrr) * cos(y*rrr + 2*pi*x*rr) ;
        % v(ii,:) = vv(xx,yy);
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

            A = chunkermat(chnkr, Dk);   
            A = A + eye(chnkr.npt);
            f = det(A);
            if nargout == 2
                derDk = 2*kernel('helm', 'fd_d', zk);  
                derA = chunkermat(chnkr, derDk);    
                varargout{1} = f*trace(derA*inv(A));
            end

        else
                
            dval = 1;
    
            F = chunkerflam(chnkr, Dk, dval, opts_flam);
            f = exp(rskelf_logdet(F));
            if nargout == 2
                A = chunkermat(chnkr, Dk);   
                A = A + eye(chnkr.npt);
                derDk = 2*kernel('helm', 'fd_d', zk);  
                derA = chunkermat(chnkr, derDk);    
                varargout{1} = f*trace(derA*inv(A));
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
                A = chunkermat(chnkr, Dk);    
                A = A + eye(chnkr.npt);
                opti = sum( ones(1, length(u))*( ( (u.').*( A\(v.') ) ).*(chnkr.wts(:)) ) );
                f = 1/opti;
                if nargout ==2 
                    Derk = 2*kernel('helmholtz','fd_d',zk);
                    derA = chunkermat(chnkr, Derk);
                    deropti = ones(1, length(u))*...
                                (   (u.').*(  A \( derA*(A\(v.')) )  ).*(chnkr.wts(:))   );
                    varargout{1} = deropti/(opti^2);
                end
            

            else
                F = chunkerflam(chnkr, Dk, 1.0, opts_flam); 
                sol = (rskelf_sv(F,v.'));
                opp = (u.').*sol.*chnkr.wts(:);
                opti = sum(opp(:));  
                f = 1/opti;
                if nargout == 2
                    Derk = 2*kernel('helmholtz','fd_d',zk);
                    choice = [];
                    choice.corrections = true;
                    cormat = chunkermat(chnkr, Derk, choice);
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
                end
            end
        end
    otherwise
        fprintf('Choice is invalid\n');
end
end