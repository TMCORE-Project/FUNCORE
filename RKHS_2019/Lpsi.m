function [dpsidlon,dpsidlat] = Lpsi(beta,dKdlon,dKdlat)
n        = size(beta,2);
dpsidlon = zeros(n,n);
dpsidlat = zeros(n,n);
for i = 1:n
    for j = 1:i
        dpsidlon(:,i) = dpsidlon(:,i) + beta(i,j) .* dKdlon(:,i);
        dpsidlat(:,i) = dpsidlat(:,i) + beta(i,j) .* dKdlat(:,i);
    end
end