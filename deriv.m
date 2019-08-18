function d = deriv(var,matrix,kdtree,a)

d = sum(matrix .* var(kdtree),2) / a;