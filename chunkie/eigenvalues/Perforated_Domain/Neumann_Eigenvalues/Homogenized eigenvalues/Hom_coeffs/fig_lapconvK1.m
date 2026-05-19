%function [I, I1] = fig_lapconvK1


% ans for m = 50, N = 250, M = 90;   is -0.0551947380820081 
% ans for m = 100, N = 500, M = 180;   is -0.0551947382407008
% ans for m = 500, N = 1000, M = 500;   is -0.0551947382394182


warning('off','MATLAB:nearlySingularMatrix')  % backward-stable ill-cond is ok!
warning('off','MATLAB:rankDeficientMatrix')
lso.RECT = true; 

format long g

jumps = [-1 0];   % jump we are after


jumps1 = [1 0];   % Alex's jump

% U.e1 = 1; U.e2 = 1i;    
% U.nei = 1; 
% [tx ty] = meshgrid(-U.nei:U.nei)
% U.trlist = tx(:)+1i*ty(:); % 3x3

U.e1 = 1; U.e2 = 1i;    
U.nei = 1; 
[tx ty] = meshgrid(-U.nei:U.nei); 
U.trlist = tx(:)+1i*ty(:); % 3x3

% m = 50; 
m = 50;
[U L R B T] = doublywalls(U,m);
proxyrep = @LapSLP;  


Rp = 1.4;
% M = 90;        
M = 90;
p.x = Rp * exp(1i*(0:M-1)'/M*2*pi); 
p = setupquad(p); 

a = 0.25; b = 0.25;  
src.x = 0.12+.13i; 
src.w = 1; 
src.nx = 1+0i; 
src.nei = 5;   

% -------------------------- single soln and plot
% N = 250; 
N = 250;
s = wormcurve(a,b,N);
rhs1 = [0*s.x; jumps1(1)+0*L.x; 0*L.x; jumps1(2)+0*B.x; 0*B.x];   % driving
rhs = [0*s.x; jumps(1)+0*L.x; 0*L.x; jumps(2)+0*B.x; 0*B.x];   % driving
tic
E = ELSmatrix(s, p, proxyrep, U);              % fill
co = linsolve(E,rhs,lso);                   % direct bkw stable solve
co1 = linsolve(E,rhs1,lso);                   % direct bkw stable solve
toc
fprintf('resid norm = %.3g\n', norm(E*co - rhs));
fprintf('resid norm1 = %.3g\n', norm(E*co1 - rhs1));
u_bdry = evalsol(s, p, proxyrep, U, s.x, co);
u_bdry1 = evalsol(s,p,proxyrep,U,s.x,co1);
xcoord = real(s.x);
ycoord = imag(s.x);
w = s.w;
I = sum(u_bdry .* xcoord .* w);
I1 = sum(u_bdry1 .* xcoord .* w);
Iy1 = sum(u_bdry .* ycoord .* w);
fprintf('Integral = %.16g\n', I);
fprintf('Integral1 = %.16g\n', I1);
fprintf('Integraly1 = %.16g\n', Iy1);

%%
sum(w) - pi/2;
intx = sum(xcoord.*w);
intxx = sum((xcoord.^2).*w);
fprintf('Integral of x*x on unit circle = %.16g and its difference from pi/64 is = %d\n', intxx, abs(intxx-pi/64));


%keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u = evalsol(s,p,proxyrep,U,z,co)

N = numel(s.x);
sig = co(1:N); psi = co(N+1:end);     
u = proxyrep(struct('x',z), p, psi);  
for i=1:numel(U.trlist)         
  u = u + srcsum(@LapSLP, U.trlist(i),[], struct('x',z), s, sig);  % pot of SLP
end
end


function [E A B C Q] = ELSmatrix(s,p,proxyrep,U)

[~,A] = srcsum(@LapSLP, U.trlist,[], s,s);   
N = numel(s.x); 
A = A - eye(N)/2;            
[~,B] = proxyrep(s,p);         
C = Cblock(s,U,@LapSLP);
[QL QLn] = proxyrep(U.L,p); 
[QR QRn] = proxyrep(U.R,p);
[QB QBn] = proxyrep(U.B,p); 
[QT QTn] = proxyrep(U.T,p);
Q = [QR-QL; QRn-QLn; QT-QB; QTn-QBn];
E = [A B; C Q];
end

function C = Cblock(s,U,densrep)     
nei = U.nei;
N = numel(s.x);
m = numel(U.L.x);
[CL CLn] = srcsum(densrep,nei*U.e1 + (-nei:nei)*U.e2,[], U.L,s);
[CR CRn] = srcsum(densrep,-nei*U.e1 + (-nei:nei)*U.e2,[], U.R,s);
[CB CBn] = srcsum(densrep,(-nei:nei)*U.e1 + nei*U.e2,[], U.B,s);
[CT CTn] = srcsum(densrep,(-nei:nei)*U.e1 - nei*U.e2,[], U.T,s);
C = [CR-CL; CRn-CLn; CT-CB; CTn-CBn];
end


function s = wormcurve(a,b,N)
% wormcurve creates the hole inside the unit cell

s.t = 2*pi*(0:N-1)'/N;
s.x = a*cos(s.t)+b*1i*sin(s.t);
% s.x = s.x + 0.3i*sin(2*real(s.x));   % spec diff can limit acc in this shape
s = setupquad(s,N);  
s.inside = @(z) (real(z)/a).^2+(imag(z)/b).^2 < 1; % yuk
end

function [U L R B T] = doublywalls(U,M)
% doublywalls creates the outer walls of the unit cell and assigns the
% normal at each wall. Observe that the normal at left wall is inward, 
% while at the right wall it's outward. Similarly,the normal at bottom 
% wall is inward, while at the top wall it's outward.


[x w] = gauss(M); 
w=w/2;
L.x = (-U.e1 + U.e2*x)/2; 
L.nx = (-1i*U.e2)/abs(U.e2) + 0*L.x; 
L.w=w*abs(U.e2);   
R = L; 
R.x = L.x + U.e1;
B.x = (-U.e2 + U.e1*x)/2; 
B.nx = (1i*U.e1)/abs(U.e1) + 0*B.x; 
B.w=w*abs(U.e1);
T = B; 
T.x = B.x + U.e2;
U.L = L; 
U.T = T; 
U.R = R; 
U.B = B;
end

function [x,w] = gauss(N)
beta = .5./sqrt(1-(2*(1:N-1)).^(-2));
T = diag(beta,1) + diag(beta,-1);
[V,D] = eig(T);
x = diag(D); 
[x,i] = sort(x);
w = 2*V(1,i).^2;
end
