ll = 1:0.1:4;
nh0 = 1;
for ii = 2:length(ll)
    [~,nh] = get_cao_chunker(ll(ii));
    if nh ~= nh0
        fprintf('ll = %d, holes = %d\n', ll(ii), nh);
        nh0 = nh;
    end
end