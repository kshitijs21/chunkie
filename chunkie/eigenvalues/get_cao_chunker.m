function [chnkr, nholes] = get_cao_chunker(j)
% this shifts the chunker from get_chunker to chunkers with reference cell
% Y centered at origin.

rad = 1; ctr = [0.0;0.0];
circfun = @(t) ctr + rad*[cos(t(:).');sin(t(:).')];

l = j;
k = 1/(4*l); % radius of hole
h = 5*l; 

rmax = 2/h;
nch = max(10, ceil(2*pi/rmax) + 2);

opts = [];
opts.maxchunklen = 1.1;
pref = []; pref.nchmax = 20000;

chnkr1 = chunkerfuncuni(circfun, nch, opts, pref);
plot(chnkr1, 'k.');
chnkreps = chnkr1;


chunkerhole0 = k*chunkerfuncuni(circfun, 6); %% hole of order 1/(4*l) centered at origin 
chunkerhole0 = chunkerhole0.reverse();       %% hole with normals poitning inwards
% if l >=1
%     nh = 1;
% end
nh = 0;
hold on,
chunkers(10000,1) = chunker();
chunkers(1) = chnkr1;
ictr = 1;
for i = -ceil(l-1) : ceil(l-1)
    for j = -ceil(l-1) : ceil(l-1)

        ctr = [i/l, j/l];
        sq = 1/2 * [1/l , 1/l];
        subsq = 1/2 * [1/l - 2/h, 1/l - 2/h];
        l1 = norm(ctr + [-1,1].*subsq);
        l2 = norm(ctr + [-1,-1].*subsq);
        l3 = norm(ctr + [1,-1].*subsq);
        l4 = norm(ctr + [1,1].*subsq);

        % % % % hor = sqrt(1-(ctr(1) - sq(1))^2); 
        % % % % ver = sqrt(1-(ctr(2) - sq(1))^2); 
        % % % % xx = -ver:0.001:ver;
        % % % % yy = -hor:0.001:hor;
        % % % % plot(ctr(1) - sq(1), yy, 'r.');
        % % % % % plot(ctr(2) - sq(1), xx, 'r.');
        % % % % plot(xx, ctr(2) - sq(1),'r.');
        % plot(xx, bb(1), 'r.');
        % % % % plot(xx, bb(2), 'r.');
        % % % plot(ctr(1) - sq(1), -1:0.01:1, 'r.');
        % % % % plot(ctr(2) - sq(1), xx, 'r.');
        % % % plot(-1:0.01:1, ctr(2) - sq(1),'r.');
        % % % 
        % % % 
        x = (ctr(1) - subsq(1)):0.01:(ctr(1) + subsq(1));
        y = (ctr(2) - subsq(1)):0.01:(ctr(2) + subsq(1));
        a = [ctr(1)-subsq(1),ctr(1)+subsq(1)];
        b = [ctr(2)-subsq(1) ctr(2)+subsq(1)];
        plot(a(1), y, 'm.');
        plot(a(2), y, 'm.');
        plot(x, b(1), 'm.');
        plot(x, b(2), 'm.');


        % check1 = [i/l,j/l] - 1/(2*l)*[1.0,1.0];               % in the first quadrant it covers the all possible vertices of epsilon-cubes that are closest to origin, assuming the unit cell has perforations centered at (1/2,1/2)
        % l1 = norm(check1 + [1/h,1/h]);    % for the epsilon-cube in the first quadrant, l1 corresponds to that vertice of sub-cube(cenetered at the center of epsilon-cube such that it's edge is of size epsilon - 2/h) which is near to check1
        % l2 = norm(check1 + [1/l,0.0] + [-1/h,1/h]);   % this is another vertices of the same sub-cube considered as in l1
        % l3 = norm(check1 + [0.0,1/l] + [1/h,-1/h]);   % this is another vertices of the same sub-cube considered as in l1
        % l4 = norm(check1 + [1/l,1/l] + [-1/h,-1/h]);  % this is another vertices of the same sub-cube considered as in l1
        % plot(check1(1) + 1/h, (check1(2) + 1/h):0.01:(check1(2) + 1/l - 1/h), 'm.');
        % plot(check1(1) + 1/l-1/h, (check1(2) + 1/h):0.01:(check1(2) + 1/l - 1/h), 'm.');
        % plot( (check1(1) + 1/h):0.01:(check1(1) + 1/l - 1/h), check1(2) + 1/h,'m.');
        % plot( (check1(1) + 1/h):0.01:(check1(1) + 1/l - 1/h), check1(2) + 1/l - 1/h,'m.');
        if (l1<=1) && (l2<=1) && (l3<=1) && (l4<=1)
            nh = nh +1;
            % plot(i/l,j/l, 'r.', 'MarkerSize', 10);

            chnkrhole = ctr + chunkerhole0;
            % chnkrctr1 = check1 + (1/(2*l))*[1.0,1.0];  % center of the epsilon cube whose vertice closer to origin is check1
            % chnkrctr2 = [-1,1].*chnkrctr1;             % this is in the quadrant of negative x and y
            % chnkrctr3 = [1,-1].*chnkrctr1;             % this is in the quadrant of x and negative y
            % chnkrctr4 = [-1,-1].*chnkrctr1;            % this is in the quadrant of negative x and negative y
            % chnkrhole1 = chnkrctr1 + chunkerhole0;     % introducing hole at chnkrctr1
            % chnkrhole2 = chnkrctr2 + chunkerhole0;     % introducing hole at chnkrctr2
            % chnkrhole3 = chnkrctr3 + chunkerhole0;     % introducing hole at chnkrctr3
            % chnkrhole4 = chnkrctr4 + chunkerhole0;     % introducing hole at chnkrctr4
            chunkers(ictr+1) = chnkrhole;
            % chunkers(ictr+2) = chnkrhole2;
            % chunkers(ictr+3) = chnkrhole3;
            % chunkers(ictr+4) = chnkrhole4;
            ictr = ictr + 1;
            % chnkreps = merge([chnkreps, chnkrhole]);
            % chnkreps = merge([chnkreps, chnkrhole1, chnkrhole2, chnkrhole3, chnkrhole4]);
            % nh = nh + 4;
           
            
        end
        % plot(check1(1), -sqrt(1-check1(1)^2):0.001:sqrt(1-check1(1)^2), 'b.');
        %     plot(-check1(1), -sqrt(1-check1(1)^2):0.001:sqrt(1-check1(1)^2), 'b.');
        % 
        %     plot(-sqrt(1-check1(2)^2):0.001:sqrt(1-check1(2)^2), check1(2), 'r.');
        %     plot(-sqrt(1-check1(2)^2):0.001:sqrt(1-check1(2)^2), -check1(2), 'r.');
    end
 end
    chnkr = merge(chunkers(1:ictr));
    % chnkr = chnkreps;
    nholes = nh;
    % hor = sqrt(1-(ctr(1) + sq(1))^2); 
    % ver = sqrt(1-(ctr(2) + sq(1))^2); 
    % xx = -ver:0.01:ver;
    % yy = -hor:0.01:hor;
    % plot(ctr(1) + sq(1), yy, 'r.');
    % plot(xx, ctr(2) + sq(1),'r.');
    % % % plot(ctr(1) + sq(1), -1:0.01:1, 'r.');
    % % % plot(-1:0.01:1, ctr(2) + sq(1),'r.');
    plot(0,0, 'r.', 'MarkerSize', 10);
    plot(chnkr, 'k.');
end
