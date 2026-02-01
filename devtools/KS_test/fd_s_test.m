%%%%%%%%%% TEST FOR fd_s %%%%%%%%%%%

clear;

chnkr = get_chunker(1);
hh = 1:1:10;
err = zeros(1,length(hh));
err2 = zeros(1,length(hh));
err1 = zeros(1,length(hh));
h0 = 1e-2;
zz = 1.82;
Freqk = kernel('helm', 'fd_s', zz);
Ak = chunkermat(chnkr, Freqk); 
for i = 1:length(hh)
     h = h0/hh(i);
     Dp = kernel('helm', 's', zz + h);  
     Ap = chunkermat(chnkr, Dp);
     Dm = kernel('helm', 's', zz - h);
     Am = chunkermat(chnkr, Dm);
     err1(i) = norm(Ap - Am);
     err(i) = norm((1/h)*(Ap - Am) - 2*Ak);
     err2(i) = norm(Ap - Am - 2*h*Ak);
end
loglog(hh, err, 'k.');
hold on; loglog(hh, (h0./hh).^2, 'b.');
% figure, loglog(hh, err2, 'm.');
% hold on; loglog(hh, (h0./hh), 'r.');
%hold on; loglog(zk, (h0./zk), 'r.');



% zk = 0:0.01:1;
% Dk = kernel('helm', 'd', 3);  
% Freqk = kernel('helm', 'freq_diff', 3);
% Ak = chunkermat(chnkr, Dk);  
% Fk = chunkermat(chnkr, Freqk);
% for i = 1:length(zk)
%     Dp = kernel('helm', 'd', 3 + zk(i));  
%     Ap = chunkermat(chnkr, Dk); 
%     Dm = kernel('helm', 'd', 3 - zk(i));
%     Am = chunkermat(chnkr, Dk);
%     Freq = kernel('helm', 'freq_diff', 3 + zk(i));
%     F = chunkermat(chnkr, Freqk);
% end
% %%
% 
%     det1 = det(A);
% end
