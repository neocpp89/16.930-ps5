function [AL, AR] = face_adjoint_consistency_tangent(left_master, right_master, left_jacobians, right_jacobians, nu)
    % the element AL(1,1) should go in A_global([left element's right node],[left element's left node]).
    % element AL(2,1) should go in A_global([left element's right node], [left element's 2nd to left node]).
    % element AR(1,1) goes in A_global([right element's left node], [right element's left node])
    w = [1; -1]; % wL - wR
    
    % hack for now...
    jl = left_jacobians(1);
    jr = right_jacobians(1);

    % this looks backwards, but we are looking at the right side of the left element, and vice versa for the right element
    AL = -left_master.dright*0.5*nu*w'./jl;
    AR = -right_master.dleft*0.5*nu*w'./jr;
end
