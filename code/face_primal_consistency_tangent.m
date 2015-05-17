function [AL, AR] = face_primal_consistency_tangent(left_master, right_master, left_jacobians, right_jacobians, nu)
    % the element AL(1,1) should go in A_global([left element's right node],[left element's left node]).
    % element AL(1,2) should go in A_global([left element's right node], [left element's 2nd to left node]).
    % element AR(1,1) goes in A_global([right element's left node], [right element's left node])
    w = [1; -1]; % wL - wR
    
    % hack for now...
    jl = left_jacobians(1);
    jr = right_jacobians(1);

    AL = -w*0.5*nu*left_master.dright'./jl;
    AR = -w*0.5*nu*right_master.dleft'./jr;
end
