function [A, f] = assemble(mesh, nu, b, c, ffn)
    master = mkmaster(mesh);

    ne = size(mesh.dgnodes, 2);
    ng = size(master.gp, 1);
    nf = size(mesh.f, 1);

    J = zeros(ng, ne);
    for i=1:ne
        J(:, i) = master.dphi'*mesh.dgnodes(:, i);
    end

    A = sparse(mesh.ndof, mesh.ndof);
    f = zeros(mesh.ndof, 1);

    % volume integrals
    for i=1:ne
        enn = mesh.nn(:, i);
        scale = master.gw;
        S = diag(scale);
        JM = diag(J(:, 1));
        IJM = diag(1 ./ J(:, 1));
        Ac = element_convection_tangent(master, J(:, 1), c);
        Anu = element_diffusion_tangent(master, J(:, 1), nu);
        Ab = element_reaction_tangent(master, J(:, 1), b);

        fg = ffn(master.phi'*mesh.dgnodes(:,i)); % evaluate f at guass points

        fk = master.phi*S*JM*fg;
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
            vv = 0.5*[c + abs(c); c - abs(c)];
            ww = [1; -1];
            nlr = [nl, nr];

            Af1 = ww*vv';
            A(nlr, nlr) = A(nlr, nlr) + Af1;

            eln = mesh.nn(:, el);
            ern = mesh.nn(:, er);

            % ugly hack to get (constant) jacobians in there.
            jl = J(1,el);
            jr = J(1,er);

            % this looks backwards, but we are looking at the right side of the left element, and vice versa for the right element
            Af2l = -ww*0.5*nu*master.dright'/jl;
            Af2r = -ww*0.5*nu*master.dleft'/jr;

            A(nlr, eln) = A(nlr, eln) + Af2l;
            A(nlr, ern) = A(nlr, ern) + Af2r;

            % adjoint consistency
            A(eln, nlr) = A(eln, nlr) + Af2l';
            A(ern, nlr) = A(ern, nlr) + Af2r';

            % same as before, looks backward but it's the right side of left element and vice versa.
            etaf = 2;
            rfl = master.ar(end)/jl;
            rfr = master.al(1)/jr;
            Af3 = -0.5*etaf*nu*(rfl+rfr)*ww*ww';
            A(nlr, nlr) = A(nlr, nlr) + Af3;
        else
            % external (boundary) face
            %{
            if (el < 0)
                % left boundary
                vv = 0.5*[c - abs(c)];
                ww = [-1];

                Af1 = ww*vv';
                % A(nr, nr) = A(nr, nr) + Af1;

                ern = mesh.nn(:, er);

                % ugly hack to get (constant) jacobians in there.
                jr = J(1,er);

                % this looks backwards, but we are looking at the right side of the left element, and vice versa for the right element
                Af2r = -ww*nu*master.dleft'/jr;

                A(nr, ern) = A(nr, ern) + Af2r;

                % adjoint consistency
                A(ern, nr) = A(ern, nr) + Af2r';

                % same as before, looks backward but it's the right side of left element and vice versa.
                etaf = 2;
                rfr = -etaf*nu*master.al/jr;
                Af3r = ww*[-rfr(1)]';
                A(nr, nr) = A(nr, nr) + Af3r;
                f(nr) = f(nr) - mesh.bcv(1)*rfr(1);
            else
                % right boundary
                vv = 0.5*[c + abs(c)];
                ww = [1];

                Af1 = ww*vv';
                % A(nl, nl) = A(nl, nl) + Af1;

                eln = mesh.nn(:, el);

                % ugly hack to get (constant) jacobians in there.
                jl = J(1,el);

                % this looks backwards, but we are looking at the right side of the left element, and vice versa for the right element
                Af2l = -ww*nu*master.dright'/jl;

                A(nl, eln) = A(nl, eln) + Af2l;

                % adjoint consistency
                A(eln, nl) = A(eln, nl) + Af2l';

                % same as before, looks backward but it's the right side of left element and vice versa.
                etaf = 2;
                rfl = -etaf*nu*master.ar/jl;
                Af3l = ww*[rfl(end)]';
                A(nl, nl) = A(nl, nl) + Af3l;
                f(nl) = f(nl) + mesh.bcv(2)*rfl(end);
            end
            %}
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

    % A = full(A);
end
