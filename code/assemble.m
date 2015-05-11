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
        scale = master.gw .* J(:, i);
        invdetJ = 1 ./ J(:, i);
        S = diag(scale);
        IDJ = diag(invdetJ);
        Ac = -master.dphi*S*IDJ*c*master.phi';
        Anu = master.dphi*S*IDJ*IDJ*nu*master.dphi';
        Ab = master.phi*S*b*master.phi';

        fg = ffn(master.phi'*mesh.dgnodes(:,i)); % evaluate f at guass points

        fk = master.phi*S*fg;
        Ak = Ac + Anu + Ab;
        f(enn) = f(enn) +  fk;
        A(enn, enn) = A(enn, enn) + Ak;
    end

    % face integrals
    for i=1:nf
        if (mesh.f(i, 4) > 0 && mesh.f(i, 3) > 0)
            % internal face
            er = mesh.f(i, 4);
            el = mesh.f(i, 3);
            nr = mesh.f(i, 2);
            nl = mesh.f(i, 1);

            Af1 = zeros(2,2);
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
            rfl = -etaf*0.5*nu*master.ar;
            rfr = -etaf*0.5*nu*master.al;
            Af3l = ww*[rfl(end); -rfl(end)]';
            Af3r = ww*[rfr(1); -rfr(1)]';
            A(nlr, nlr) = A(nlr, nlr) + Af3l + Af3r;

        else
            % external (boundary) face
            if (mesh.f(i, 3) < 0)
                % left boundary
            else
                % right boundary
            end
        end
    end

    % boundaries
    lambdaL = mesh.bcnn(1);
    wL = mesh.lrn(1);
    A(wL, lambdaL) = 1;
    A(lambdaL, wL) = 1;
    lambdaR = mesh.bcnn(2);
    wR = mesh.lrn(2);
    A(wR, lambdaR) = (-1);
    A(lambdaR, wR) = 1;
    f(mesh.bcnn) = mesh.bcv;

    % A = full(A);
end
