[chnkr, nh] = get_chunker(2);
xx = chnkr.r(1,:);
yy = chnkr.r(2,:);

u = (xx.^2 - yy.^2 + yy).*sin(xx) + cos(xx).*(0.5*yy.^2 + 0.3*yy);
v = (xx - yy.^3 + xx.^2).*cos(-yy) + sin(yy).*(xx.^2 + xx);


z0 = 3.409274696865717;       
n = 10;
for i = 1:n
    [zfast, yfast, zslow, yslow] = update_iterates(z0, chnkr, u, v);
    z0 = zfast;
    z0slow = zslow;
    fprintf('i = %d, Difference in eigenvalues using fast and slow code = %d\n', i, abs(zfast-zslow));
    if (abs(yfast) < 1e-12)
        break;
    end
end
niters = i;
fprintf(' Number of holes = %d\n Iterations required = %d\n Eignvalue = %d\n Optimisation function value = %d', nh, i, zfast, yfast)



function [z, y, zslow, yslow] = update_iterates(zk, chnkr1, u, v)
    chnkr = chnkr1;
    Dkslow = 2*kernel('helm', 'd', zk);  
    Aslow = chunkermat(chnkr1, Dkslow);    
    Aslow = Aslow + eye(chnkr1.npt);
    derDkslow = 2*kernel('helm', 'freq_diff', zk);  
    derAslow = chunkermat(chnkr1, derDkslow);
    optislow = ones(1, length(u))*( ( (u.').*( Aslow\(v.') ) ).*(chnkr.wts(:)) );
    deroptislow = ones(1, length(u))*...
                (   (u.').*(  Aslow \( derAslow*(Aslow\(v.')) )  ).*(chnkr.wts(:))   );
    yslow = 1/optislow;
    zslow = zk - optislow/deroptislow;

    
    Dk = 2*kernel('helm', 'd', zk); 
    F = chunkerflam(chnkr1, Dk, 1.0);
    derDk = 2*kernel('helm', 'freq_diff', zk);  
    derA = chunkermat(chnkr1, derDk);  
    opti = ones(1, length(u))*( ( (u.').*( rskelf_sv(F,v.') ) ).*(chnkr1.wts(:)) );
    deropti = ones(1, length(u))*...
                  (   (u.').*(  rskelf_sv(F, ( derA*(rskelf_sv(F,v.')) ) )  ).*(chnkr1.wts(:))   );
    y = 1/opti;
    z = zk - opti/deropti;
end