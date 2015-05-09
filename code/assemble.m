function [A, f] = assemble(mesh, nu, b, c, ffn)
    master = mkmaster(mesh);

    ne = size(mesh.dgnodes, 2);
    ng = size(master.gp, 1);
    nf = size(mesh.f, 1);

    J = zeros(ng, ne);
    for i=1:ne
        J(:, i) = master.dphi'*mesh.dgnodes(:, i)
    end

    A = sparse(mesh.ndof, mesh.ndof);
    f = zeros(mesh.ndof, 1);

    % volume integrals
    for i=1:ne
        enn = mesh.nn(:, i);
        scale = master.gw .* J(:, i);
        S = diag(scale);
        Ac = -master.dphi*S*c*master.phi';
        Anu = master.dphi*S*nu*master.dphi';
        Ab = master.phi*S*b*master.phi';

        fg = ffn(master.phi'*mesh.dgnodes(:,i)); % evaluate f at guass points

        fk = master.phi*S*fg;
        Ak = Ac + Anu + Ab;
        f(enn) = f(enn) +  fk;
        A(enn, enn) = A(enn, enn) + Ak;
    end

    % face integrals
    for i=1:nf
        if (mesh.f(i, 4) > 0)
            % internal face
        else
            % external (boundary) face
        end
    end

    % boundaries
    % for i=1:numel(mesh.bcnn)
    % end

    A = full(A);
end
