chnkr = get_chunker(3);
xx = chnkr.r(1,:); yy = chnkr.r(2,:);

u1 = (xx.^2 - yy.^2 + yy).*sin(xx) + cos(xx).*(0.5*yy.^2 + 0.3*yy);
v1 = (xx - yy.^3 + xx.^2).*cos(-yy) + sin(yy).*(xx.^2 + xx);
u2 = (xx.^2 - yy.^2 + sin(yy)).*sin(xx) + cos(xx).*(0.5*sin(yy).^2 + 0.3*yy);
v2 = (sin(xx) - yy.^3 + xx.^2).*cos(-yy) + sin(yy).*(xx.^2 + sin(xx));
u3 = (xx.^2 + yy).*sin(xx) + cos(xx).*(0.3*yy);
v3 = (xx + xx.^2).*cos(-yy) + sin(yy).*(xx.^2 + xx);

z0 = 3.409274696865717;  
n = 10;                
f = @(zk) update_iterates(zk, chnkr, u1, v1, u2, v2, u3, v3);
for i = 1:n
    [z, y] = f(z0);
    z0 = z;
    fprintf('i=%d,  f=%d,  z0=%d\n',i, abs(y), real(z));
    if (abs(y) < 1e-12)
        break;
    end
end
niters = i;


function [z, y] = update_iterates(zk, chnkr1, u1, v1, u2, v2, u3, v3)

    Dk = 2*kernel('helm', 'd', zk); 
    F = chunkerflam(chnkr1, Dk, 1.0);
    derDk = 2*kernel('helm', 'fd_d', zk);  
    derA = chunkermat(chnkr1, derDk);  
    opti = ones(1, length(u1))*( ( (u1.').*( rskelf_sv(F,v1.') ) ).*(chnkr1.wts(:)) ) ...
             + ones(1, length(u2))*( ( (u2.').*( rskelf_sv(F,v2.') ) ).*(chnkr1.wts(:)) )...
             + ones(1, length(u1))*( ( (u3.').*( rskelf_sv(F,v3.') ) ).*(chnkr1.wts(:)) );
    deropti = ones(1, length(u1))*...
                  (   (u1.').*(  rskelf_sv(F, ( derA*(rskelf_sv(F,v1.')) ) )  ).*(chnkr1.wts(:))   )...
              + ones(1, length(u2))*...
                  (   (u2.').*(  rskelf_sv(F, ( derA*(rskelf_sv(F,v2.')) ) )  ).*(chnkr1.wts(:))   )...
              + ones(1, length(u3))*...
                  (   (u3.').*(  rskelf_sv(F, ( derA*(rskelf_sv(F,v3.')) ) )  ).*(chnkr1.wts(:))   )   ;
    y = 1/opti;
    z = zk - opti/deropti;
end