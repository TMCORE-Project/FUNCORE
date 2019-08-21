function [w,idx] = gen_weights_2d (pts_i,pts_c,eps,poly_deg,nSamples,base_opt)
% Input parameters
% pts_i Column vectors with stencil node locations; approximation to
% be accurate at pts_c
% eps Power of r in RBF fi(r) = r^m, with m odd, >= 3.
% poly_deg Degree of supplementary polynomials (d = -1 no polynomials)
% base_opt, 1 for r^m base, 2 for Gaussian base
% Output parameter
% w Matrix with three columns, containing weights for d/dlon, d/dlat

x  = pts_i(:,1);
y  = pts_i(:,2);
xc = pts_c(:,1);
yc = pts_c(:,2);

x  = x - xc;
y  = y - yc;

dist = sqrt(x.^2 + y.^2);

[~,idx_tmp] = sort(dist);
rd          = dist(idx_tmp(1:nSamples));
idx         = idx_tmp(1:nSamples);

xd = x(idx); % Chosen points x coordinate
yd = y(idx); % Chosen points x coordinate

pts_d = pts_i(idx,:);

n = length(x);%x

% ------ RBF part --------------------------------------------------------
r_1d = pdist(pts_d);
r    = squareform(r_1d);

% Choose base function
if base_opt == 1
    phi      = r.^eps;
    dphidx   = -eps .* rd.^(eps-2) .* xd;
    dphidy   = -eps .* rd.^(eps-2) .* yd;
elseif base_opt == 2
    phi      = exp(-(eps.*r).^2);
    dphidx   = 2 * eps^2.*exp(-(eps.* rd).^2) .* xd;
    dphidy   = 2 * eps^2.*exp(-(eps.* rd).^2) .* yd;
end

L0 = [dphidx,dphidy]; % RHSs

% ------ Polynomial part -------------------------------------------------
if poly_deg == -1 % Special case; no polynomial terms,
    A = phi;
    L = L0; % i.e. pure RBF
else % Create matrix with polynomial terms and matching constraints
    X = xd(:,ones(1,poly_deg+1)); X(:,1) = 1; X = cumprod( X,2);
    Y = yd(:,ones(1,poly_deg+1)); Y(:,1) = 1; Y = cumprod( Y,2);
    
    np       = (poly_deg+1)*(poly_deg+2)/2; % Number of polynomial terms
    XY       = zeros(nSamples,np);
    col      = 1; % Assemble polynomial matrix block
    for k = 0:poly_deg
        XY(:,col:col+k) = X(:,k+1:-1:1).*Y(:,1:k+1);
        col = col+k+1;
    end
    
    L1 = zeros(np,2); % Create matching RHSs
    if poly_deg >= 1; L1(2,1) = 1; L1(3,2) = 1; end
    
    A = [phi,XY;XY',zeros(np)]; % Assemble linear system to be solved
    L = [L0;L1]; % Assemble RHSs
end

% ------ Solve for weights -----------------------------------------------
W = A\L;
w = W(1:nSamples,:); % Extract the RBF-FD weights

% save mtx_my
