function [c] = loading_saved_data(z0, nh0, x, y)
% Input arguments are:
%   nh0 = hole for which you know the data corresponding to some ll less than a
%   z0 = initial guess corresponding to nh0
%   a = the value from where the data you have to extract
%   b = the value till where the data you have to extract
% Output arguments are:
%   c = struct consisting all the data in the saved file so far

ll = x:0.01:y;
zz = z0;
nhole = nh0;
counter = 0;
data(length(ll)).zero0 = []; 
data(length(ll)).fnit = [];
data(length(ll)).it = [];
data(length(ll)).nit = [];
data(length(ll)).val0 = [];
data(length(ll)).timeit = []; 
data(length(ll)).nholes = []; 
for ii = 1:length(ll)
    [~, nh] = get_chunker(ll(ii));
    if nh~=nhole
        fprintf('Saving data for Hole = %d, ll = %d', nh, ll(ii));
        counter = counter + 1;
        filename = ['ir_with_' num2str(nh) '_holes_and_initial_guess_' num2str(zz) '.mat'];
        a = load(filename);
        zz = real(a.zero0);
        nhole = nh;
        data(counter).zero0 = a.zero0;
        data(counter).fnit = a.fnit;
        data(counter).it = a.it;
        data(counter).nit = a.nit;
        data(counter).val0 = a.val0;
        data(counter).timeit = a.timeit;
        data(counter).nholes = nh;
    end 
end
c = data(1:counter);
end