% Defn of chunker
clear;
chnkr = get_chunker(3);
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
    u(ii,:) = (xx.^(2*ii + 2) - yy.^(2*ii + 2) + yy.^(2*ii + 3)).*sin(xx) + cos(xx).*(0.5*yy.^(2*ii + 2) + 0.3*yy);
    v(ii,:) = (xx.^(2*ii + 3) - yy.^(3*ii + 2) + cos(xx.^(2*ii + 2))).*cos(-yy) + sin(yy).*(xx.^(2*ii + 2) + xx);
end

u(nn+1,:) = (xx.^2 - yy.^2 + yy).*sin(xx) + cos(xx).*(0.5*yy.^2 + 0.3*yy);
v(nn+1,:) = (xx - yy.^3 + xx.^2).*cos(-yy) + sin(yy).*(xx.^2 + xx);

% Initial Guess

z0 = 3.409274686896810;  

% Newton 
n = 10;              % number of iterations
f = @(zk) update_iterates(zk, chnkr, u, v, u(nn+1,:), v(nn+1,:));
for i = 1:n
    [z, y, z2, y2] = f(z0);
    z0 = z;
    fprintf('i=%d, fm=%d, fs=%d, zm=%d, zs=%d, diffz=%d, diffy=%d\n', i, abs(y2), abs(y), real(z2), real(z), abs(z-z2), abs(y-y2));
    if (abs(y2) < 1e-15)
        break;
    end
end
niters = i;

% Newton iterates calculation

function [z, y, z2, y2] = update_iterates(zk, chnkr, u, v, u1, v1)
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

    % [rr,~] = size(u);
    % opti = 0;
    % deropti = 0;
    % for ii = 1:rr
    %     sol = (rskelf_sv(F,v(ii,:).'));
    %     u_eval = chunkerkerneval(chnkr,dkern,sol,chnkr,opts);
    %     opp = ones(1, length(u(ii,:)))*( ( (u(ii,:).').*( sol ) ).*(chnkr.wts(:)) );
    %     opti = opti + opp;
    %     deropp = ones(1, length(u(ii,:)))*(   (u(ii,:).').*(  rskelf_sv(F, u_eval)  ).*(chnkr.wts(:))   );
    %     deropti = deropti + deropp;
    % end

    %%%%%%%%  SOME ELEGANT CODING  %%%%%%%%%

    % for multi optimisation function
    
    sol = (rskelf_sv(F,v.'));
    u_eval = chunkerkerneval(chnkr,dkern,sol,chnkr,opts);
    opp = u*(sol.*chnkr.wts(:));
    solder = rskelf_sv(F, u_eval);
    deropp = u*(solder.*chnkr.wts(:));
    opti = trace(opp);
    deropti = trace(deropp);
    y2 = 1/opti;
    z2 = zk - opti/deropti;

    % for single optimisation function

    sol1 = (rskelf_sv(F,v1.'));
    u_eval1 = chunkerkerneval(chnkr,dkern,sol1,chnkr,opts);
    opp1 = u1*(sol1.*chnkr.wts(:));
    solder1 = rskelf_sv(F, u_eval1);
    deropp1 = u1*(solder1.*chnkr.wts(:));
    y = 1/opp1;
    z = zk - opp1/deropp1;
end