%
% This is a playground for testing RBF and RKHS function expansion.
%

clc
clear

test_case_1d = 'cardinal';
num_node_1d = 21;

test_case_2d = 'const_one';
num_node_2d = 100;

% 1D interpolation
if strcmp(test_case_1d, 'noisy_sine')
  xi = linspace(0, 2 * pi, num_node_1d);
  yi = sin(xi)' + 0.2 * rand(num_node_1d, 1);
  plot_y_range = [-1.5 1.5];
elseif strcmp(test_case_1d, 'cardinal')
  xi = linspace(-1, 1, num_node_1d);
  yi = zeros(num_node_1d, 1);
  yi(mean(1:size(xi, 2))) = 1;
  plot_y_range = [-0.5 1.5];
elseif strcmp(test_case_1d, 'const_one')
  xi = linspace(-1, 1, num_node_1d);
  yi = ones(num_node_1d, 1);
  plot_y_range = [0.995 1.005];
end

fig = figure();

subplot(4, 2, 1)
test_rbf_1d(xi, yi, plot_y_range, 'kernel', 'ga', 'eps', 3, 'poly_deg', -1)

subplot(4, 2, 2)
test_rbf_1d(xi, yi, plot_y_range, 'kernel', 'ga', 'eps', 3, 'poly_deg', 2)

subplot(4, 2, 3)
test_rbf_1d(xi, yi, plot_y_range, 'kernel', 'ga', 'eps', 3, 'poly_deg', 1)

subplot(4, 2, 4)
test_rbf_1d(xi, yi, plot_y_range, 'kernel', 'ga', 'eps', 9, 'poly_deg', -1)

subplot(4, 2, 5)
test_rbf_1d(xi, yi, plot_y_range, 'kernel', 'ga', 'eps', 11, 'poly_deg', -1)

subplot(4, 2, 6)
test_rbf_1d(xi, yi, plot_y_range, 'kernel', 'ga', 'eps', 11, 'poly_deg', 0)

subplot(4, 2, 7)
test_rbf_1d(xi, yi, plot_y_range, 'kernel', 'phs', 'pow', 3, 'poly_deg', -1)

subplot(4, 2, 8)
test_rbf_1d(xi, yi, plot_y_range, 'kernel', 'phs', 'pow', 3, 'poly_deg', 1)

% -------------------------------------------------------------------------

function res = dist(x1, x2)

  res = norm(x2 - x1);

end

function res = kernel_ga(x, x0, varargin)

  eps = cell2mat(varargin(4));
  res = exp(-(dist(x, x0) * eps)^2);

end

function res = kernel_phs(x, x0, varargin)

  pow = cell2mat(varargin(4));
  res = (dist(x, x0))^pow;

end

function [indices, col] = create_polynomial(d, l, col, indices)

  if ~exist('indices', 'var')
    indices = zeros(l, nchoosek(l + d, l) - 1);
    col = 1;
    for i = 1 : l
      indices(1:i,col) = 1;
      col = col + 1;
      [indices, col] = create_polynomial(d, i, col, indices);
    end
  else
    while 1
      j = col - 1;
      done = 1;
      for i = l : -1 : 1
        if indices(i,j) ~= d
          done = 0;
          break
        end
      end
      if done
        return
      end
      indices(:,col) = indices(:,j);
      indices(i:l,col) = indices(i,col) + 1;
      col = col + 1;
    end
  end

end

function res = create_coef_matrix(x, kernel)

  n = size(x, 2);
  res = zeros(n, n); % Note: Coefficient matrix is full.
  for i = 1 : n
    for j = 1 : n
      res(i,j) = kernel(x(i), x(j));
    end
  end

end

function P = create_polynomial_matrix(x, l)

  [d, n] = size(x);
  pi = create_polynomial(d, l);

  P = ones(n, size(pi, 2) + 1);
  for i = 1 : n
    for j = 2 : size(P, 2)
      for k = 1 : d
        if pi(k,j-1) ~= 0
          P(i,j) = P(i,j) * x(pi(k,j-1),i);
        end
      end
    end
  end
  
end

function res = interp(xi, kernel, w, xo)

  n = size(xi, 2);
  res = 0;
  for i = 1 : n
    res = res + w(i) * kernel(xi(i), xo);
  end

end

function test_rbf_1d(xi, yi, plot_y_range, varargin)

  kernel_name = cell2mat(varargin(2));
  display(['==> Using kernel_', kernel_name])
  if strcmp(kernel_name, 'ga')
    kernel = @(x, x0) kernel_ga(x, x0, varargin{:});
  elseif strcmp(kernel_name, 'phs')
    kernel = @(x, x0) kernel_phs(x, x0, varargin{:});
  end

  [d, n] = size(xi);
  A = create_coef_matrix(xi, kernel);
  l = cell2mat(varargin(6));
  xo = linspace(min(xi) - 0.2 * (max(xi) - min(xi)), max(xi) + 0.2 * (max(xi) - min(xi)), 201);

  if l == -1
    w = A \ yi;
    yo = arrayfun(@(x) interp(xi, kernel, w, x), xo);
  else
    P = create_polynomial_matrix(xi, l);
    A = [A P; P' zeros(size(P, 2))];
    L = [yi; zeros(nchoosek(l + d, l), 1)];
    W = A \ L;
    w = W(1:n);
    yo = arrayfun(@(x) interp(xi, kernel, w, x), xo);
    P = create_polynomial_matrix(xo, l);
    yo = yo + W(n+1:end)' * P';
  end

  plot(xo, yo, 'r-', 'linewidth', 3);
  hold on;
  plot(xi, yi, 'o', 'markerfacecolor', 'blue');

  set(gca, 'linewidth', 4, 'fontsize', 16)
  xlim([min(xo) max(xo)])
  ylim(plot_y_range)
  grid on
  if strcmp(kernel_name, 'ga')
    eps = cell2mat(varargin(4));
    if l == -1
      title(['Gaussian RBF with shape parameter ', num2str(eps)])
    else
      title(['Gaussian RBF with polynomials and shape parameter ', num2str(eps)])
    end
  elseif strcmp(kernel_name, 'phs')
    if l == -1
      title('PHS RBF without polynomials')
    else
      title('PHS RBF with polynomials')
    end
  end

end
