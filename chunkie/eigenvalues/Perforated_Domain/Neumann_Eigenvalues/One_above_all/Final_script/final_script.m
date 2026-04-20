% ll = 5:0.01:15;
ll = [2 2.6 3];
opts = [];
[~,nh0] = get_chunker(ll(1));
% z0 = 3.449551391360910;  % initial guess corresponding to ll = 5 or number of hole = 68
z0 = 3.449763720345;
func = 'ir';
method = 's';
fprintf('Objective function = %s\n', func);
fprintf('Method = %s\n', method);

opts.newtiter = 15;
for ii = 2:length(ll)
    [chnkr, nh] = get_chunker(ll(ii));
    opts.newtguess = z0;
    opts.secguess = [z0 - 0.05, z0];
    if nh ~= nh0
        fprintf('Hole = %d: Initial Guess = %d\n', nh, z0);
        [zero0, val0, nit, valit, timeit] = zeros_fred_opti(func, method, chnkr, opts);
        fprintf('  Root = %d, fn value = %d\n', zero0, val0);
        it = valit(1,:);
        fnit = valit(2,:);
        save(['Secant/' char(func) '_with_' num2str(nh) '_holes_and_initial_guess_' num2str(z0) '.mat'], 'zero0', 'val0', 'nit', 'it', "fnit", "timeit"); 
        z0 = real(zero0);
    end
    nh0 = nh;
end