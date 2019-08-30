function [e,beta] = orthogonalization(a)
% Stabilized Gram Schmidt method for orthonormalization
n    = size(a,2);
b    = zeros(n,n);
e    = zeros(n,n);
beta = zeros(n,n);
for i = 1:n
   b(:,i) = a(:,i);
   for j = 1:i-1
       beta(i,j) = dot(b(:,j),a(:,i)) / dot(b(:,j),b(:,j));
       b   (:,i) = b(:,i) - beta(i,j) * b(:,j);
   end
   e(:,i) = b(:,i) / norm(b(:,i));
end
