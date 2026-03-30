function derf = der_fred_opti(zk, chnkr, choice, obj)
% fred_opti returns either the Fredholm determinant and its derivative or
% the optimisation function and its derivative based on the value of 'choice'.
%
% 'choice' is a character such that it takes the following values
%       * 'F' for Fredholm determinant and it's derivative
%       * 'O' for Optimisation function (1/<u_i,inv(A)v_i>) and its derivative
% 
% 'obj' is a struct with the entries 'u' and 'v' such that 
%       * obj.u is a (nd, chnkr.npt) matrix
%       * obj.v is a (nd, chnkr.npt) matrix
%       * nd is provided by the user

if strcmp(choice, 'F')

    Dk = 2*kernel('helm', 'd', zk);  
    A = chunkermat(chnkr, Dk);   
    A = A + eye(chnkr.npt);
    f = det(A);
    derDk = 2*kernel('helm', 'fd_d', zk);  
    derA = chunkermat(chnkr, derDk);    
    derf = f*trace(derA*inv(A));

elseif strcmp(choice, 'O')
    u = obj.u;
    v = obj.v;
    % u = input('Enter the vector u:\n');
    % v = input('Enter the vector v:\n');
    [ru, cu] = size(u);
    [rv, cv] = size(v);
    if cu~=chnkr.npt || cv~=chnkr.npt || cu~=cv || ru~=rv
        fprintf('Dimensions of u and v are not correct\n');
    else
    Dk = 2*kernel('helm', 'd', zk); 
    F = chunkerflam(chnkr, Dk, 1.0);
    dkern = 2*kernel('helmholtz','fd_d',zk);
    choice = [];
    choice.corrections = true;
    cormat = chunkermat(chnkr, dkern, choice);
    choice = [];
    choice.forcesmooth = true;
    choice.cormat = cormat;

    sol = (rskelf_sv(F,v.'));
    u_eval = zeros(chnkr.npt, ru);
    for ii = 1:ru
        u_eval(:,ii) = chunkerkerneval(chnkr,dkern,sol(:,ii),chnkr,choice);
    end

    opp = (u.').*sol.*chnkr.wts(:);            %%% MATRIX of u_i A^{-1} v_j    
    solder = rskelf_sv(F, u_eval);             %%% MATRIX of u_i (A^{-1})' v_j
    deropp = (u.').*(solder.*chnkr.wts(:));
    opti = sum(opp(:));                        %%% sum of u_i A^{-1} v_i     
    deropti = sum(deropp(:));
    f = 1/opti;
    derf = deropti/(opti^2);
    end
else
    fprintf('choice is invalid\n');
end
end