%%%%%%%%%%% eigval computation for neumann laplacian using the complex
%%%%%%%%%%% muller for Fred determinant
clear;
m = 2;

% f1 = @(j,zk) get_determinant(zk, get_chunker(j));
x = zeros(m,1);
y = zeros(m,1); 
niters = zeros(m,1);
ll = 2;

for j = 1:length(ll)
    [chnkr, nh] = get_chunker(ll(j));
    f = @(zk) get_determinant(zk, chnkr);

    if j==1
        x0 = 3.21; 
        x1 = 3.24;
        x2 = 3.27;
    else
        x0 = x(j-1) - 0.01;
        x1 = x(j-1);
        x2 = x(j-1) + 0.01;
    end
    y0 = f(x0); y1 = f(x1); y2 = f(x2);     % 3 function evaluation for each j
    t1 = tic;
    if (abs(y0) < 1e-12)
        x(j) = x0;
    elseif (abs(y1) < 1e-12)
        x(j) = x1;
    elseif (abs(y2) < 1e-12)
        x(j) = x2;
    else 
        nit = 0;
        while(abs(y2) > 1e-12)
            nit = nit+1;
            h0 = x1 - x0; 
            h1 = x2 - x1;
            delta0 = (y1 - y0)/h0;
            delta1 = (y2 - y1)/h1;
            a = (delta1 - delta0)/(h1 + h0); 
            b = a*h1 + delta1;
            c = y2;

            x3 = x2 - 2*c/(b + sign(b)*sqrt(b*b - 4*a*c));

            x0 = x1; y0 = y1;
            x1 = x2; y1 = y2;
            x2 = x3; y2 = f(x2);   % 1 function evaluation for each i 
        end
        niters(j) = nit;
        x(j) = x2;
        y(j) = y2;
    end
    time = toc(t1);
    fprintf('Holes = %d, Time = %d, Iterations = %d, fn eval = %d, root = %d.\n', nh, time, niters(j), 3 + niters(j), x(j));
end

function det1 = get_determinant(zk, chnkr1)

    Dk = 2*kernel('helm', 'd', zk);  
    A = chunkermat(chnkr1, Dk);   
%%
    A = A + eye(chnkr1.npt);
%%
    det1 = det(A);
end
