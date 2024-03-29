function w = gen_weights (x,y,z,m,d)
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
coord    = [x,y,z];
r_1d     = pdist(coord);
r        = squareform(r_1d);
% phi      = r.^m; % RBF matrix
% dphidr   = m .* r(:,1).^(m-2); % dphidr = dphidr / r to avoide divde 0;
phi      = exp( - m^2 * r.^2 );
dphidx   = -2 * m^2 .* phi(:,1) .* xd;
dphidy   = -2 * m^2 .* phi(:,1) .* yd;
dphidz   = -2 * m^2 .* phi(:,1) .* zd;
L0       = [dphidx,dphidy,dphidz]; % RHSs
% ------ Polynomial part -------------------------------------------------
if d == -1 % Special case; no polynomial terms,
    A = phi;
    L = L0; % i.e. pure RBF
else % Create matrix with polynomial terms and matching constraints
    X = xd(:,ones(1,d+1)); X(:,1) = 1; X = cumprod( X,2);
    Y = yd(:,ones(1,d+1)); Y(:,1) = 1; Y = cumprod( Y,2);
    Z = zd(:,ones(1,d+1)); Z(:,1) = 1; Z = cumprod( Z,2);
    
    np       = (d+1)*(d+2)*(d+3)/6; % Number of polynomial terms
    XYZ      = zeros(n,np);
    XYZ(:,1) = 1;
    for k = 1:d
        ids  =  k   *(k+1)*(k+2)/6+1;
        ide  = (k+1)*(k+2)*(k+3)/6;
        XYZ(:,ids:ids+k) = X(:,k+1:-1:1) .* Y(:,1:k+1);
        
        idz1 = (k-1)*(k  )*(k+1)/6+1;
        idz2 = (k  )*(k+1)*(k+2)/6;
        XYZ(:,ids+k+1:ide) = XYZ(:,idz1:idz2) .* Z(:,2);
    end
    
    L1 = zeros(np,2); % Create matching RHSs
    if d >= 1
        L1(2,1) = 1;
        L1(3,2) = 1;
        L1(4,3) = 1;
    end
    A = [phi,XYZ;XYZ',zeros(np)]; % Assemble linear system to be solved
    L = [L0;L1]; % Assemble RHSs
end
% ------ Solve for weights -----------------------------------------------
W = A\L;
w = W(1:n,:); % Extract the RBF-FD weights
