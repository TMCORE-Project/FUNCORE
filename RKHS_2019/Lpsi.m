function [dpsidlon,dpsidlat] = Lpsi(beta,dKdlon,dKdlat)
n        = size(beta,2);
dpsidlon = zeros(n,1);
dpsidlat = zeros(n,1);
for i = 1:n
    for j = 1:i
        dpsidlon(i,1) = dpsidlon(i,1) + beta(i,j) .* dKdlon(i,1);
        dpsidlat(i,1) = dpsidlat(i,1) + beta(i,j) .* dKdlat(i,1);
    end
end