function [P] = lagrange(xpts)
    % code modified from 2006-11-20 Dan Ellis dpwe@ee.columbia.edu
    N = length(xpts);
    P = zeros(N,N);
    for i = 1:N
      pp = poly(xpts((1:N) ~= i));
      P(i,:) = pp ./ polyval(pp, xpts(i));
    end
end
