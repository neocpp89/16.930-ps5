function [f] = element_load_vector(master, jacobians, f_gauss)
    % f_gauss has the load function evaluated at the gauss points for this element.
    S = diag(master.gw .* jacobians);
    f = master.phi*S*f_gauss;
end
