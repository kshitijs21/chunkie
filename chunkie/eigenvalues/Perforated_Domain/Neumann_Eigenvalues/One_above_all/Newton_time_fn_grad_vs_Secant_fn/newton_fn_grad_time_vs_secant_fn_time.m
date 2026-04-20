l = 3;
[chnkr, nh] = get_chunker(l);
func = 'ir';
method = 'n';
z00 = 3.396302810422;
opts = [];
opts.newtguess = z00;
opts.secguess = [z00-0.1 z00];
fprintf('Newton: \n');
[rt, frt, nit, valit, timen] = zeros_fred_opti(func, method, chnkr, opts);
timenewt = timen(1,:) + timen(2,:);

fprintf('Secant: \n');
mthd = 's';
[rts, frts, nits, valits, times] = zeros_fred_opti(func, mthd, chnkr, opts);
m = min(length(timenewt), length(times));
timediff = timenewt(1:1:m) - times(1:1:m);