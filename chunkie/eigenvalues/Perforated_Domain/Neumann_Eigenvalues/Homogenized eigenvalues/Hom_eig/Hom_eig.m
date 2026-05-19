clear;
chnkr = get_chunker(1);
func = 'fred';
method = 'c';
opts = [];
opts.cheb1 = 2.8;
opts.cheb2 = 4;
[root, val, f, rts] = zeros_fred_opti(func, method, chnkr, opts);