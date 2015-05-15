function [A] = element_diffusion_tangent(master, jacobians, b)
    % multiply by J for integration.
    S = diag(master.gw .* jacobians);
    A = master.phi*S*b*master.phi';
end
