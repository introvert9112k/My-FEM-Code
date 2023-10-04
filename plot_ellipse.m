function plot_ellipse(X, varargin)
% Plots an ellipse of the form x'*X*x <= 1

% plot vectors of the form a*u + b*v where u,v are eigenvectors of X
X = [0.8 0;0 0.4];
[V,D] = eig(X);
u = V(:,1);
v = V(:,2);
l1 = D(1,1);
l2 = D(2,2);

pts = [];

delta = .01;

for alpha = -1/sqrt(l1)-delta:delta:1/sqrt(l1)+delta
    beta = sqrt((1 - alpha^2 * l1)/l2);
    pts(:,end+1) = alpha*u + beta*v;
end
for alpha = 1/sqrt(l1)+delta:-delta:-1/sqrt(l1)-delta
    beta = -sqrt((1 - alpha^2 * l1)/l2);
    pts(:,end+1) = alpha*u + beta*v;
end

plot(pts(1,:), pts(2,:), varargin{:})
end