function w = gen_weights_2d (x,y,m,d,base_opt)
% Input parameters
% x,y,z Column vectors with stencil node locations; approximation to
% be accurate at x(1),y(1),z(1)
% m Power of r in RBF fi(r) = r^m, with m odd, >= 3.
% d Degree of supplementary polynomials (d = -1 no polynomials)
% base_opt, 1 for r^m base, 2 for Gaussian base
% Output parameter
% w Matrix with three columns, containing weights for d/dlon, d/dlat

xd = x-x(1);
yd = y-y(1);% Shift nodes so stencil centered at origin

n = length(xd);%x

% ------ RBF part --------------------------------------------------------
coord    = [x,y];
r_1d     = pdist(coord);
r        = squareform(r_1d);

% Choose base function
if base_opt == 1
    phi      = r.^m; % RBF matrix
    dphidx   = -m .* r(1,:)'.^(m-2) .* xd;
    dphidy   = -m .* r(1,:)'.^(m-2) .* yd;
elseif base_opt == 2
    phi      = exp(-(m.*r).^2);
    dphidx   = -2 * m^2.*exp(-(m.* r(1,:)').^2) .* xd;
    dphidy   = -2 * m^2.*exp(-(m.* r(1,:)').^2) .* yd;
end

L0 = [dphidx,dphidy]; % RHSs

% ------ Polynomial part -------------------------------------------------
if d == -1 % Special case; no polynomial terms,
    A = phi;
    L = L0; % i.e. pure RBF
else % Create matrix with polynomial terms and matching constraints
    X = xd(:,ones(1,d+1)); X(:,1) = 1; X = cumprod( X,2);
    Y = yd(:,ones(1,d+1)); Y(:,1) = 1; Y = cumprod( Y,2);
    
    np       = (d+1)*(d+2)/2; % Number of polynomial terms
    XY       = zeros(n,np);
    col      = 1; % Assemble polynomial matrix block
    for k = 0:d
        XY(:,col:col+k) = X(:,k+1:-1:1).*Y(:,1:k+1);
        col = col+k+1;
    end
    
    L1 = zeros(np,2); % Create matching RHSs
    if d >= 1; L1(2,1) = 1; L1(3,2) = 1; end
    
    %     if d >= 1
    %         L1(2,1) = 1;
    %         L1(3,2) = 1;
    %     end
    %     if d >= 2
    %         L1(4 ,1) = 2 * x(1);
    %         L1(5 ,1) = y(1);
    %         L1(5 ,2) = x(1);
    %         L1(6 ,2) = 2 * y(1);
    %     end
    
    A = [phi,XY;XY',zeros(np)]; % Assemble linear system to be solved
    L = [L0;L1]; % Assemble RHSs
end

% ------ Solve for weights -----------------------------------------------
W = A\L;
w = W(1:n,:); % Extract the RBF-FD weights

% save mtx_my
