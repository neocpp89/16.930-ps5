function [A] = face_dg_flux_tangent(master, jacobians, c)
    % the element A(1,1) should go in A_global([left element's right node],[left element's right node]).
    % element A(2,2) should go in A_global([right element's left node], [right element's left node]).
    % A(1,2) and A(2,1) follow the same pattern
    v = 0.5 * [c + abs(c); c - abs(c)];
    w = [1; -1]; % wL - wR
    A = w*v';
end
