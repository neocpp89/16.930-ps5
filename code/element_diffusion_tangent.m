function [A] = element_diffusion_tangent(master, jacobians, nu)
    % multiply by J for integration, divide by J^2 for derivatives of functions.
    S = diag(master.gw ./ jacobians);
    A = master.dphi*S*nu*master.dphi';
end
