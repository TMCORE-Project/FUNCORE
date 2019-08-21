function w = RBF_FD_PHS_pol_weights (x,y,m,d)
% Input parameters
% x,y Column vectors with stencil node locations; approximation to
% be accurate at x(1),y(1)
% m Power of r in RBF fi(r) = r^m, with m odd, >= 3.
% d Degree of supplementary polynomials (d = -1 no polynomials)
%
% Output parameter
% w Matrix with three columns, containing weights for d/dx, d/dy,
% and the Laplacian d2/dx2+d2/dy2, respectively.
x = x-x(1); y = y-y(1); % Shift nodes so stencil centered at origin
n = length(x);%x
% ------ RBF part --------------------------------------------------------
A0 = hypot(bsxfun(@minus,x,x'),bsxfun(@minus,y,y')).^m; % RBF matrix
L0 = m*(bsxfun(@times,(hypot(x,y)).^(m-2),[-x,-y,m*ones(n,1)])); % RHSs
% ------ Polynomial part -------------------------------------------------
if d == -1 % Special case; no polynomial terms,
A = A0; L = L0; % i.e. pure RBF
else % Create matrix with polynomial terms and matching constraints
X = x(:,ones(1,d+1)); X(:,1) = 1; X = cumprod( X,2);
Y = y(:,ones(1,d+1)); Y(:,1) = 1; Y = cumprod( Y,2);
np = (d+1)*(d+2)/2; % Number of polynomial terms
XY = zeros(n,np); col = 1; % Assemble polynomial matrix block
for k = 0:d
XY(:,col:col+k) = X(:,k+1:-1:1).*Y(:,1:k+1);
col = col+k+1;
end
L1 = zeros(np,3); % Create matching RHSs
if d >= 1; L1(2,1) = 1; L1(3,2) = 1; end
if d >= 2; L1(4,3) = 2; L1(6,3) = 2; end
A = [A0,XY;XY',zeros(col-1)]; % Assemble linear system to be solved
L = [L0;L1]; % Assemble RHSs
end
% ------ Solve for weights -----------------------------------------------
W = A\L;
w = W(1:n,:); % Extract the RBF-FD weights

% save mtx_old