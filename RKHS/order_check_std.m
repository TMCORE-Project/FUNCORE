clc
clear

eps       = 1;
nSamples  = 7;

res(1) = 0.25;
res(2) = 0.5;
res(3) = 1;
res(4) = 2;
res(5) = 4;
res(6) = 8;
res(7) = 16;

d2r = pi/180;
r2d = 180/pi;

res = res*d2r;

for ires = 1:size(res,2)
    lon = 0:res(ires):2*pi-res(ires);
    lon = lon';
    lat = zeros(size(lon));
    
    np  = size(lon,1);
    
    dist   = distance(0.,0.,lat,lon,'radians');
    dist2d = bsxfun(@circshift,dist,0:size(dist,1)-1);
    
    xd              = dist;
    xd(end/2+1:end) = -xd(end/2+1:end);
    xd              = bsxfun(@circshift,xd,0:size(xd,1)-1)';
    
    [v_tmp,idx_tmp] = sort(dist2d(1,:));
    v       = xd(idx_tmp(1:nSamples));
    idx     = idx_tmp(1:nSamples);
    weights = rbfga_weights('x',eps,v',dist(1));
    
    s0 = sin(v);
    c0 = cos(v);
    d0 = sum(s0.*weights');
    
    L1(ires) = sum(abs(d0-c0(1))) / sum(abs(c0(1)));
    L2(ires) = sqrt( sum((d0-c0(1)).^2)./sum(c0(1).^2) );
    
    if ires >= 2
        order_L1(ires-1) = log(L1(ires)/L1(ires-1))/log(res(ires)/res(ires-1));
        order_L2(ires-1) = log(L2(ires)/L2(ires-1))/log(res(ires)/res(ires-1));
    end
end
