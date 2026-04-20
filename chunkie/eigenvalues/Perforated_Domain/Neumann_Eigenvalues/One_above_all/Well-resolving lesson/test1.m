ll = 2;
[chnkr, nh] = get_chunker(ll);
chnkr = sort(chnkr);
xx = chnkr.r(1,:);
yy = chnkr.r(2,:);
func = 'ir';
method = 'c';

% z00 = 3.235377405;
% a = 3.235377405 - 0.2;
% b = 3.235377405 + 0.2;
% opts = [];
% opts.nd = 4;
% opts.cheb1 = a;
% opts.cheb2 = b;
% 
% opts.newtguess = 2.8;
% opts.secguess = [2.8 2.9];
% opts.cmguess = [2.7 2.8 2.9];
% xz = fred_opti(zk, chnkr, func, opts);
% [zero0, val, nit, rts] = zeros_fred_opti(func, method, chnkr, opts);
% 
% return

% rr = rand(length(xx),1);
% rrr = rand(length(xx),1);
% %nd = 1;
% ii = 5;
% uu = @(x, y) rr + x*rrr + (x.*(y.^ii))*rr + (x.^2)*rrr + (y.^ii)*rr) * sin(2*pi*x*rr) ;
% u = uu(xx,yy);
% u = zeros(nd, length(xx));
% v = zeros(nd, length(xx));
% uu = zeros(nd, length(xx));
% vv = zeros(nd, length(xx));
%rr = rand(6,1);
%rr2 = -1 + 2*rand(3,1);

% rr = -1+2*rand(6,1);
% rr2 = -1 + 2*rand(3,1);
% 
 rr = 1 + rand(6,1);
        
        uu = @(x, y) (rr(1)*x + rr(2)*y + rr(3)*(x.*y) + rr(4)*x.^2 + rr(5)*y.^2 + rr(6)) .* sin(pi*x*rr(1)) ;
        
        for ii = 2:2
            u = uu(xx,yy);
        end

%%%% 
% rr = rand(6,1);
% rr2 = rand(3,1);
% 
% u = @(x,y) (rr(1) + rr(2)*x + rr(3)*y + rr(4)*x.^2 + rr(5)*x.*y + rr(6).*y.^2).*sin(rr2(1) + rr2(2)*x + rr2(3)*y);
 % for ii = 2:2
 %        uu(1,:) = (xx.^(2*(1/(ii)) + 2) - yy.^(2*(1/(ii)) + 2) + yy.^(2*(1/(ii)) + 3)).*sin(xx) + cos(xx).*(0.5*yy.^(2*(1/(ii)) + 2) + 0.3*yy);
 %        vv(1,:) = (xx.^(2*(1/(ii)) + 3) - yy.^(3*(1/(ii)) + 2) + cos(xx.^(2*(1/(ii)) + 2))).*cos(-yy) + sin(yy).*(xx.^(2*(1/(ii)) + 2) + xx);
 % end
 % 



u = u(:);
%v = v(:);

%uu = uu(:);
%vv = vv(:);

% well resolved test for u


figure(2)
clf
plot(u, 'k.')



ww = chnkr.wts(:);
figure(3)
clf
plot(cumsum(ww(:)), u, 'k.');

[~, ~, umat, ~] = lege.exps(chnkr.k);

uuse = reshape(u, 16, []);
ucoefs = umat*uuse;

uerr = sqrt(ucoefs(15,:).^2 + ucoefs(16,:).^2);
figure(4)
clf
semilogy(abs(uerr), 'k.');

return
figure(5)
clf
plot(uu, 'k.')

www = chnkr.wts(:);
figure(6)
clf
plot(cumsum(www(:)), uu, 'k.');

[~, ~, uumat, ~] = lege.exps(chnkr.k);

uuuse = reshape(uu, 16, []);
uucoefs = uumat*uuuse;

uuerr = sqrt(uucoefs(15,:).^2 + uucoefs(16,:).^2);
figure(7)
clf
semilogy(abs(uuerr), 'k.');


% return

% well resolved test for v


figure(2)
clf
plot(v, 'k.')



ww = chnkr.wts(:);
figure(3)
clf
plot(cumsum(ww(:)), v, 'k.');

[~, ~, vmat, ~] = lege.exps(chnkr.k);

vuse = reshape(u, 16, []);
vcoefs = vmat*vuse;

verr = sqrt(vcoefs(15,:).^2 + vcoefs(16,:).^2);
figure(4)
clf
semilogy(abs(verr), 'k.');


% figure(5)
% clf
% plot(uu, 'k.')
% 
% www = chnkr.wts(:);
% figure(6)
% clf
% plot(cumsum(www(:)), uu, 'k.');
% 
% [~, ~, uumat, ~] = lege.exps(chnkr.k);
% 
% uuuse = reshape(uu, 16, []);
% uucoefs = uumat*uuuse;

% uuerr = sqrt(uucoefs(15,:).^2 + uucoefs(16,:).^2);
% figure(7)
% clf
% semilogy(abs(uuerr), 'k.');
