function [A, f] = assemble(mesh, nu, b, c, ffn)
    master = mkmaster(mesh);

    ne = size(mesh.dgnodes, 2);
    ng = size(master.gp, 1);
    nf = size(mesh.f, 1);

    J = master.dphi'*mesh.dgnodes;

    A = sparse(mesh.ndof, mesh.ndof);
    f = zeros(mesh.ndof, 1);

    fg = ffn(master.phi'*mesh.dgnodes); % evaluate f at guass points

    % volume integrals
    for i=1:ne
        enn = mesh.nn(:, i);

        Ac = element_convection_tangent(master, J(:, i), c);
        Anu = element_diffusion_tangent(master, J(:, i), nu);
        Ab = element_reaction_tangent(master, J(:, i), b);
        fk = element_load_vector(master, J(:, i), fg(:, i));

        Ak = Ac + Anu + Ab;
        f(enn) = f(enn) +  fk;
        A(enn, enn) = A(enn, enn) + Ak;
    end

    % face integrals
    for i=1:nf
        nl = mesh.f(i, 1);
        nr = mesh.f(i, 2);
        el = mesh.f(i, 3);
        er = mesh.f(i, 4);
        if (el  > 0 && er > 0)
            % internal face
            ww = [1; -1];
            nlr = [nl, nr];
            elr = [el, er];

            Af1 = face_dg_flux_tangent(master, J(:, elr), c);
            A(nlr, nlr) = A(nlr, nlr) + Af1;

            eln = mesh.nn(:, el);
            ern = mesh.nn(:, er);

            [AL, AR] = face_primal_consistency_tangent(master, master, J(:, el), J(:, er), nu);
            A(nlr, eln) = A(nlr, eln) + AL;
            A(nlr, ern) = A(nlr, ern) + AR;

            [AL, AR] = face_adjoint_consistency_tangent(master, master, J(:, el), J(:, er), nu);
            A(eln, nlr) = A(eln, nlr) + AL;
            A(ern, nlr) = A(ern, nlr) + AR;

            % same as before, looks backward but it's the right side of left element and vice versa.
            Af3 = face_lifting_tangent(master, master, J(:, el), J(:, er), nu);
            A(nlr, nlr) = A(nlr, nlr) + Af3;
        else
            % external (boundary) face (do this later)
        end
    end

    % boundaries
    lambdaL = mesh.bcnn(1);
    wL = mesh.lrn(1);
    A(wL, lambdaL) = -1;
    A(lambdaL, wL) = -1;
    lambdaR = mesh.bcnn(2);
    wR = mesh.lrn(2);
    A(wR, lambdaR) = 1;
    A(lambdaR, wR) = 1;
    fbc = mesh.bcv;
    fbc(1) = -fbc(1);
    f(mesh.bcnn) = fbc;
end
