function [mesh] = mkmesh_uniform(endpts, nelements, porder)
    mesh.porder = porder;
    mesh.plocal = linspace(0, 1, porder+1)'; % use uniform spacing inside elements too
    xx = linspace(endpts(1), endpts(2), nelements+1);
    h_element = xx(2) - xx(1);
    for n=1:nelements
        mesh.dgnodes(:, n) = xx(n) + h_element*mesh.plocal;
    end

    mesh.nn = zeros(size(mesh.dgnodes)); % node number index mapping
    mesh.nn(:) = 1:numel(mesh.nn);
    mesh.bcnn = [1 2] + numel(mesh.nn); % node numbers for bcs
    mesh.bcv = [0 1]; % values for BC (you can override this outside of mkmesh)
    mesh.lrn = [1 numel(mesh.nn)]; % left and right node numbers

    mesh.f = zeros(nelements+1, 4);
    mesh.f(1:nelements-1, 1) = mesh.nn(end, 1:end-1);
    mesh.f(1:nelements-1, 2) = mesh.nn(1, 2:end);
    mesh.f(1:nelements-1, 3) = 1:nelements-1;
    mesh.f(1:nelements-1, 4) = 2:nelements;
    mesh.f(end-1, :) = [mesh.bcnn(1), 1, -1, 1];
    mesh.f(end, :) = [mesh.nn(end,end), mesh.bcnn(2), nelements, -2];

    mesh.ndof = numel(mesh.nn) + 2;
end
