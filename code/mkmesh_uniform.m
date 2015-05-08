function [mesh] = mkmesh_uniform(endpts, nelements, porder)
    mesh.porder = porder;
    mesh.plocal = linspace(0, 1, porder+1)'; % use uniform spacing inside elements too
    xx = linspace(endpts(1), endpts(2), nelements+1);
    xx = xx(1:end-1);
    h_element = xx(2) - xx(1);
    for n=1:nelements
        mesh.dgnodes(:, n) = xx(n) + h_element*mesh.plocal;
    end
end
