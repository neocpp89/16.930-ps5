function [x, w] = golub_welsch(n)
    % taken from http://www.mathworks.com/matlabcentral/fileexchange/23972-chebfun/content/chebfun/examples/quad/html/GaussQuad.html
    beta = .5./sqrt(1-(2*(1:n-1)).^(-2)); % 3-term recurrence coeffs
    T = diag(beta,1) + diag(beta,-1);     % Jacobi matrix
    [V,D] = eig(T);                       % Eigenvalue decomposition
    x(:,1) = diag(D); [x,i] = sort(x);         % Legendre points
    w(:,1) = 2*V(1,i).^2;                      % Quadrature weights
end
