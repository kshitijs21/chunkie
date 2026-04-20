% rad = 1; ctr = [0.0;0.0];
% circfun = @(t) ctr + rad*[cos(t(:).');sin(t(:).')];
% 
% nch = 20;
% % nch = max(10, ceil(2*pi/rmax) + 2);
% optsc = [];
% optsc.maxchunklen = 0.02;
% prefc = []; 
% prefc.nchmax = 10000;
% chnkr = chunkerfuncuni(circfun, nch, optsc, prefc);
ll = 1:0.5:2;
zk = 3;
for ii = 1:length(ll)

    [chnkr, nh] = get_chunker(ll(ii));
    fprintf('Hole = %d \n', nh);

    % derDk = 2*kernel('helm', 'fd_d', zk);  
    %             derA = chunkermat(chnkr, derDk);    

    t1 = tic;
    derDk = 2*kernel('helm', 'fd_d', zk); 
    timee = toc(t1);
    t2 = tic;
    derA = chunkermat(chnkr, derDk);  
    timeks = toc(t2);
    % tic, cc = det(derA); toc;
    fprintf('kernel build up time = %d, chunkermat time = %d\n', timee, timeks);
end