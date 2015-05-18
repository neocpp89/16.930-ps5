function [A] = face_lifting_tangent(left_master, right_master, left_jacobians, right_jacobians, nu)
    % the element A(1,1) should go in A_global([left element's right node],[left element's right node]).
    % element A(2,2) should go in A_global([right element's left node], [right element's left node]).
    % A(1,2) and A(2,1) follow the same pattern
    w = [1; -1]; % wL - wR
    v = [1; -1]; % vL - vR
    
    % hack for now...
    jl = left_jacobians(1);
    jr = right_jacobians(1);

    % this looks backwards, but we are looking at the right side of the left element, and vice versa for the right element
    etaf = 2;
    rfl = left_master.ar(end)/jl;
    rfr = right_master.al(1)/jr;
    A = -0.5*etaf*nu*(rfl+rfr)*w*v';
end
