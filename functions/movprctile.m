function outvec = movprctile(input, window, percentile)
%movprctile calculates moving percentile
%   outvec = movprctile(input, window, percentile)

% Output
outvec = zeros(size(input));
l = length(outvec);

for i = 1 : l
    if i < window
        outvec(i) = prctile(input(1:i), percentile);
    else
        outvec(i) = prctile(input(i-window+1:i), percentile);
    end
end

end

