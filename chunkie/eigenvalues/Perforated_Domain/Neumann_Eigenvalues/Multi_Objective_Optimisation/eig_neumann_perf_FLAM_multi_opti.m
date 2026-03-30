% Defn of chunker
clear;
chnkr = get_chunker(2.6);
xx = chnkr.r(1,:);
yy = chnkr.r(2,:);

% Defn of u and v

nn = 50;
u = zeros(nn+1, length(xx));
v = zeros(nn+1, length(xx));
% ii = 1:1:n+1;
% u = (xx.^(2*ii + 2) - sin(yy.^(2*ii + 2)) + cos(yy.^(2*ii + 3))).*sin(xx) + cos(xx).*(0.5*yy.^(2*ii + 2) + 0.3*yy);
% v = (xx.^(2*ii + 3) - exp(yy.^(3*ii + 2)) + cos(xx.^(2*ii + 2))).*cos(-yy) + sin(yy).*(xx.^(2*ii + 2) + xx);

for ii = 1:nn
    u(ii,:) = (xx.^(2*(1/(ii+1)) + 2) - yy.^(2*(1/(ii+1)) + 2) + yy.^(2*(1/(ii+1)) + 3)).*sin(xx) + cos(xx).*(0.5*yy.^(2*(1/(ii+1)) + 2) + 0.3*yy);
    v(ii,:) = (xx.^(2*(1/(ii+1)) + 3) - yy.^(3*(1/(ii+1)) + 2) + cos(xx.^(2*(1/(ii+1)) + 2))).*cos(-yy) + sin(yy).*(xx.^(2*(1/(ii+1)) + 2) + xx);
end

u(nn+1,:) = (xx.^2 - yy.^2 + yy).*sin(xx) + cos(xx).*(0.5*yy.^2 + 0.3*yy);
v(nn+1,:) = (xx - yy.^3 + xx.^2).*cos(-yy) + sin(yy).*(xx.^2 + xx);

func = 'ir';
method = 'n';
opts = [];
opts.newtguess = 3.409274686896810;
opts.newtiter = 10;
opts.u = u;
opts.v = v;

[zk, y0, it, zzz] = zeros_fred_opti(func, method, chnkr, opts);
% % Initial Guess
% return
z0 = 3.409274686896810;  

% Newton 
n = 10;              % number of iterations
zzo = zeros(1,n);
yyo = zeros(1,n);
f = @(zk) update_iterates(zk, chnkr, u, v, u(nn+1,:), v(nn+1,:));
for i = 1:n
    % fprintf('For i=%d iteration\n', i);
    [z, y, z1, y1, z2, y2] = f(z0);
    z0 = z;
    zzo(i) = z2;
    yyo(i) = y2;
    %fprintf('Elegant fme=%d, zme=%d, fse=%d, zse=%d, diffz=%d, diffy=%d\n', abs(y2), real(z2), abs(y), real(z), abs(z-z2), abs(y-y2));
    fprintf('i = %d, Elegant fme=%d, zme=%d\n', i, abs(y2), real(z2));
    % fprintf('Diff btw Elegant and Non-Elegant diffzme=%d, diffyme=%d\n', abs(z1-z2), abs(y1-y2));
    if (abs(y2) < 1e-12)
        break;
    end
end
niters = i;

% Newton iterates calculation

function [z, y, z1, y1, z2, y2] = update_iterates(zk, chnkr, u, v, u1, v1)
    Dk = 2*kernel('helm', 'd', zk); 
    F = chunkerflam(chnkr, Dk, 1.0);
    dkern = 2*kernel('helmholtz','fd_d',zk);
    opts = [];
    opts.corrections = true;
    cormat = chunkermat(chnkr, dkern, opts);
    opts = [];
    opts.forcesmooth = true;
    opts.cormat = cormat;
    %%%%%%%%  SOME NON-ELEGANT CODING  %%%%%%%%%

    [rr,~] = size(u);
    opti = 0;
    deropti = 0;
    
    for ii = 1:rr
        sol = (rskelf_sv(F,v(ii,:).'));
        u_eval = chunkerkerneval(chnkr,dkern,sol,chnkr,opts);
        opp = ones(1, length(u(ii,:)))*( ( (u(ii,:).').*( sol ) ).*(chnkr.wts(:)) );
        opti = opti + opp;
        deropp = ones(1, length(u(ii,:)))*(   (u(ii,:).').*(  rskelf_sv(F, u_eval)  ).*(chnkr.wts(:))   );
        deropti = deropti + deropp;
    end
    y1 = 1/opti;
    z1 = zk - opti/deropti;
    opti0 = opti;
    deropti0 = deropti;

%%%%%%    ELEGANT CODING  %%%%%%%%

% for multi objective function

    sol = (rskelf_sv(F,v.'));
    for ii = 1:rr
        u_eval(:,ii) = chunkerkerneval(chnkr,dkern,sol(:,ii),chnkr,opts);
    end

    opp = u*(sol.*chnkr.wts(:));    %%% MATRIX of u_i A^{-1} v_j    
    solder = rskelf_sv(F, u_eval);  %%% MATRIX of u_i (A^{-1})' v_j
    deropp = u*(solder.*chnkr.wts(:));
    opti1 = trace(opp);              %%% sum of u_i A^{-1} v_i     
    deropti1 = trace(deropp);        %%% sum of u_i (A^{-1})' v_i
    oppm = (u.').*sol.*chnkr.wts(:);    %%% MATRIX of u_i A^{-1} v_j    
    solderm = rskelf_sv(F, u_eval);  %%% MATRIX of u_i (A^{-1})' v_j
    deroppm = (u.').*(solderm.*chnkr.wts(:));
    opti1m = sum(oppm(:));              %%% sum of u_i A^{-1} v_i     
    deropti1m = sum(deroppm(:));
    % diffoptim = abs(opti1-opti1m);
    % diffderoptim = abs(deropti1-deropti1m);
    y2 = 1/opti1m;
    z2 = zk - opti1m/deropti1m;
    %fprintf('Diffopti=%d, Diffderopti=%d\n', abs(opti1m-opti0)/abs(opti0), abs(deropti1m-deropti0)/abs(deropti0));

    % for single objective function

    sol1 = (rskelf_sv(F,v1.'));
    u_eval1 = chunkerkerneval(chnkr,dkern,sol1,chnkr,opts);
    opp1 = u1*(sol1.*chnkr.wts(:));
    solder1 = rskelf_sv(F, u_eval1);
    deropp1 = u1*(solder1.*chnkr.wts(:));
    y = 1/opp1;
    z = zk - opp1/deropp1;
end