function [chnkr, nholes] = get_chunker(j)
rad = 1; ctr = [0.0;0.0];
circfun = @(t) ctr + rad*[cos(t(:).');sin(t(:).')];



l = j;
k = 1/(4*l);
h = 5*l; 

rmax = 2/h;
nch = max(10, ceil(2*pi/rmax) + 2);

opts = [];
opts.maxchunklen = 1.1;
pref = []; pref.nchmax = 20000;

chnkr1 = chunkerfuncuni(circfun, nch, opts, pref);
plot(chnkr1, 'k.');
chnkreps = chnkr1;




% opts = [];
% opts.maxchunklen = 1.1;
% pref = []; 
% pref.nchmax = 100000;
% chnkr1 = chunkerfunc(circfun, opts, pref);
% 
% l = j;
% chnkreps = chnkr1;
% nh = 0;
% k = 1/(4 * l); 
% h = 5 * l;
chunkerhole0 = k*chunkerfuncuni(circfun, 6);   %% hole of order 1/(4*l) 
chunkerhole0 = chunkerhole0.reverse();
nh = 0;
chunkers(10000,1) = chunker();
chunkers(1) = chnkr1;
ictr = 1;
hold on,
for i = 0 : l-1
    fprintf('i = %d:\n', i);
    for j = 0 : l-1
        
        fprintf('  j = %d\n', j);
        check1 = [i/l,j/l];
        sq = [1/h,1/h];
        x = (check1(1) + sq(1)):0.01:(check1(1) + 1/l - sq(1));
        y = (check1(2) + sq(1)):0.01:(check1(2) + 1/l - sq(1));
        a = [check1(1)+ sq(1),check1(1) + 1/l - sq(1)];
        b = [check1(2) + sq(1),check1(2) + 1/l - sq(1)];
        plot(a(1), y, 'm.');
        plot(a(2), y, 'm.');
        plot(x, b(1), 'm.');
        plot(x, b(2), 'm.');
        l1 = norm(check1 + [1/h,1/h]); 
        l2 = norm(check1 + [1/l,0.0] + [-1/h,1/h]);
        l3 = norm(check1 + [0.0,1/l] + [1/h,-1/h]);
        l4 = norm(check1 + [1/l,1/l] + [-1/h,-1/h]);
        if (l1<=1) && (l2<=1) && (l3<=1) && (l4<=1)
            chnkrctr1 = check1 + (1/(2*l))*[1.0,1.0];
            chnkrctr2 = [-1,1].*chnkrctr1;
            chnkrctr3 = [1,-1].*chnkrctr1;
            chnkrctr4 = [-1,-1].*chnkrctr1;
            chnkrhole1 = chnkrctr1 + chunkerhole0;
            chnkrhole2 = chnkrctr2 + chunkerhole0;
            chnkrhole3 = chnkrctr3 + chunkerhole0;
            chnkrhole4 = chnkrctr4 + chunkerhole0;
            % chunkers(ictr+1) = chnkrhole1;
            % chunkers(ictr+2) = chnkrhole2;
            % chunkers(ictr+3) = chnkrhole3;
            % chunkers(ictr+4) = chnkrhole4;
            % ictr = ictr + 4;
            chnkreps = merge([chnkreps, chnkrhole1, chnkrhole2, chnkrhole3, chnkrhole4]);
            nh = nh + 4;
         end
    end
 end
    chnkr = merge(chunkers(1:ictr));
    chnkr = chnkreps;
    nholes = nh;
    
    plot(chnkr, 'k.');
end
