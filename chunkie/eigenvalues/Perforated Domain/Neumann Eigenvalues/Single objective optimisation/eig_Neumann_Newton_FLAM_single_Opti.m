chnkr = get_chunker(2);
xx = chnkr.r(1,:);
yy = chnkr.r(2,:);

u = (xx.^2 - yy.^2 + yy).*sin(xx) + cos(xx).*(0.5*yy.^2 + 0.3*yy);
v = (xx - yy.^3 + xx.^2).*cos(-yy) + sin(yy).*(xx.^2 + xx);

z0 = 3.409274696865717;  
n = 10;% number of iterations
f = @(zk) update_iterates(zk, chnkr, u, v);
for i = 1:n
    [z, y, z2, y2] = f(z0);
    z0 = z;
    fprintf('i=%d,  f=%d,  z0=%d  diffz=%d  diffy =%d\n',i, abs(y2), real(z), abs(z-z2), abs(y-y2));
    if (abs(y) < 1e-12)
        break;
    end
end
niters = i;


function [z, y, z2, y2] = update_iterates(zk, chnkr1, u, v)
    chnkr = chnkr1;
    Dk = 2*kernel('helm', 'd', zk); 
    F = chunkerflam(chnkr1, Dk, 1.0);
    derDk = 2*kernel('helm', 'fd_d', zk);  
    derA = chunkermat(chnkr1, derDk);  
    sol = (rskelf_sv(F,v.'));
    opti = ones(1, length(u))*( ( (u.').*( sol ) ).*(chnkr1.wts(:)) );
    deropti = ones(1, length(u))*...
                  (   (u.').*(  rskelf_sv(F, ( derA*sol ) )  ).*(chnkr1.wts(:))   );
    y = 1/opti;
    z = zk - opti/deropti;

    dkern = 2*kernel('helmholtz','fd_d',zk);
    opts = [];
    opts.corrections = true;
    cormat = chunkermat(chnkr1, dkern, opts);
    opts = [];
    opts.forcesmooth = true;
    opts.cormat = cormat;
    u_eval_cor = chunkerkerneval(chnkr1,dkern,sol,chnkr1,opts);
    deropti = ones(1, length(u))*(   (u.').*(  rskelf_sv(F, u_eval_cor)  ).*(chnkr1.wts(:))   );
    y2 = 1/opti;
    z2 = zk - opti/deropti;

end


% % % % function [z, y] = update_iterates(zk, chnkr1, u, v)
% % % %     chnkr = chnkr1;
% % % %     % Dk = 2*kernel('helm', 'd', zk); 
% % % %     % F = chunkerflam(chnkr1, Dk, 1.0);
% % % %     % opti = ones(1, length(u))*( ( (u.').*( rskelf_sv(F,v.') ) ).*(chnkr1.wts(:)) );
% % % %     % sol = rskelf_sv(F,v.');
% % % %     % dkern = kernel('helmholtz','fd_d',zk);
% % % %     % opts = [];
% % % %     % opts.corrections = true;
% % % %     % cormat = chunkermat(chnkr, dkern, opts);
% % % %     % % cormat = chunkerkernevalmat(chnkr,dkern,chnkr,opts);
% % % %     % opts = [];
% % % %     % opts.forcesmooth = true;
% % % %     % opts.cormat = cormat;
% % % %     % u_eval_cor = chunkerkerneval(chnkr1,dkern,sol,chnkr1,opts);
% % % %     % deropti = ones(1, length(u))*(   (u.').*(  rskelf_sv(F, u_eval_cor)  ).*(chnkr1.wts(:))   );
% % % %     % y = 1/opti;
% % % %     % z = zk - opti/deropti;
% % % % 
% % % % end




% function [F] = get_matrix(zk, chnkr1)
% 
%     Dk = 2*kernel('helm', 'd', zk);  
%     % A = chunkermat(chnkr1, Dk);
%     F = chunkerflam(chnkr1, Dk, 1.0);
% % %%
% %     A = A + eye(chnkr1.npt);
% end
% 
% 
% function [derA] = get_dermatrix(zk, chnkr1)
%     derDk = 2*kernel('helm', 'fd_d', zk);  
%     derA = chunkermat(chnkr1, derDk);   
% end
