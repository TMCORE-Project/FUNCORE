clc
clear

eps       = 0.6;
base_type = 1;
poly_deg  = 4;

res(1) = 0.25;
res(2) = 0.5;
res(3) = 1;
res(4) = 2;
res(5) = 4;

% Select base function
if base_type == 1
    rbf_base = @(r,eps) exp(-(eps.*r).^2);
    rbf_base_deriv = @(r,eps,xd) -2*eps^2.*exp(-(eps.*r).^2).*xd;
elseif base_type == 2
    rbf_base = @(r,eps) r.^eps;
    rbf_base_deriv = @(r,eps,xd) eps * r.^(eps-2).*xd;
elseif base_type == 3
    rbf_base = @(r,eps) sqrt(1+(eps.*r).^2);
    rbf_base_deriv = @(r,eps,xd) eps.^2 ./ sqrt(1+(eps*r).^2).*xd;
end

d2r = pi/180;
r2d = 180/pi;

res = res*d2r;

for ires = 1:size(res,2)
    lon = 0:res(ires):2*pi-res(ires);
    lon = lon';
    lat = zeros(size(lon));
    
    np  = size(lon,1);
    
    f0  = sin(lon);
    df0 = cos(lon);
    
    dist   = distance(0.,0.,lat,lon,'radians');
    dist2d = bsxfun(@circshift,dist,0:size(dist,1)-1);
    
    A0  = rbf_base(dist2d,eps);
    
    xd              = dist;
    xd(end/2+1:end) = -xd(end/2+1:end);
    xd              = bsxfun(@circshift,xd,0:size(xd,1)-1)';
    
    B0 = rbf_base_deriv(dist2d,eps,xd);
    
    % Add polynominal
    if poly_deg ~= -1
        poly      = zeros(np,poly_deg);
        poly(:,1) = 1;
        for ipoly = 2:poly_deg+1
            poly(:,ipoly) = dist.^(ipoly-1);
        end
        
        A = [A0,poly;poly',zeros(poly_deg+1)];
        
        Lp = zeros(poly_deg+1,np);
        if(poly_deg>=1)
            Lp(2,:) = 1;
        end
        if(poly_deg>=2)
            for ip = 3:poly_deg
                Lp(ip,:) = ip*lon.^(ip-1);
            end
        end
        
        B = [B0;Lp];
    else
        A = A0;
        B = B0;
    end
    
    Dp = A\B;
    D  = Dp(1:np,:);
    
    df = D'*f0;
    
    dd(ires) = df(end/2+1);
    
    L1(ires) = sum(abs(df-df0)) / sum(abs(df));
    L2(ires) = sqrt( sum((df-df0).^2)./sum(df0.^2) );
    
    if ires >= 2
        order_L1(ires-1) = log(L1(ires)/L1(ires-1))/log(res(ires)/res(ires-1));
        order_L2(ires-1) = log(L2(ires)/L2(ires-1))/log(res(ires)/res(ires-1));
    end
end

order = log(abs(dd(2)+1)/abs(dd(1)+1))/log(2)