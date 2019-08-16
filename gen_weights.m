function w = gen_weights (x,y,z,lon,lat,radius,m,d)
% Input parameters
% x,y,z Column vectors with stencil node locations; approximation to
% be accurate at x(1),y(1),z(1)
% m Power of r in RBF fi(r) = r^m, with m odd, >= 3.
% d Degree of supplementary polynomials (d = -1 no polynomials)
%
% Output parameter
% w Matrix with three columns, containing weights for d/dlon, d/dlat

xd = x-x(1);
yd = y-y(1);
zd = z-z(1);% Shift nodes so stencil centered at origin

n = length(xd);%x

% ------ RBF part --------------------------------------------------------
coord = [x,y,z];
r_1d  = pdist(coord);
r     = squareform(r_1d);
A0    = r.^m; % RBF matrix
Llon  = cos(lat) * cos(lat(1)) .* sin(lon-lon(1)) * m .* r(:,1).^(m-2);
Llat  = ( sin(lat) * cos(lat(1)) .* cos(lon-lon(1)) - cos(lat) * sin(lat(1)) ) * m .* r(:,1).^(m-2);
L0    = [Llon,Llat]; % RHSs
% ------ Polynomial part -------------------------------------------------
if d == -1 % Special case; no polynomial terms,
    A = A0; L = L0; % i.e. pure RBF
else % Create matrix with polynomial terms and matching constraints
    X = xd(:,ones(1,d+1)); X(:,1) = 1; X = cumprod( X,2);
    Y = yd(:,ones(1,d+1)); Y(:,1) = 1; Y = cumprod( Y,2);
    Z = zd(:,ones(1,d+1)); Z(:,1) = 1; Z = cumprod( Z,2);
    
    np  = (d+1)*(d+2)/2; % Number of polynomial terms
    XYZ = zeros(n,np);
    col = 1; % Assemble polynomial matrix block
    for i = 0:d
        for k = 0:d
            XYZ(:,col:col+k) = X(:,k+1:-1:1) .* Y(:,1:k+1) .* Z(:,i+1);
            col = col+k+1;
        end
    end
    L1 = zeros(np,2); % Create matching RHSs
    if d >= 1
        dxdlon = - radius * sin(lon(1)) * cos(lat(1));
        dydlon = radius * cos(lon(1)) * cos(lat(1));
        dzdlon = 0;
        dxdlat = - radius * cos(lon(1)) * sin(lat(1));
        dydlat = - radius * sin(lon(1)) * sin(lat(1));
        dzdlat = radius * cos(lat(1));
        L1(2,1) = dxdlon;
        L1(2,2) = dxdlat;
        L1(3,1) = dydlon;
        L1(3,2) = dydlat;
        L1(4,1) = dzdlon;
        L1(4,2) = dzdlat;
    end
    A = [A0,XYZ;XYZ',zeros(col-1)]; % Assemble linear system to be solved
    L = [L0;L1]; % Assemble RHSs
end
% ------ Solve for weights -----------------------------------------------
W = A\L;
w = W(1:n,:); % Extract the RBF-FD weights