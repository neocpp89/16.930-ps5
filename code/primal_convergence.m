clear all; close all;

N = 2.^[4:7]';
p = [1, 2, 3]';

params = { {1, 0, 0, @(x) x.*x, 'Poisson'} };
exact_solution = {@(x) x.*(13 - x.^3)/12};
errors = zeros(numel(params), numel(p), numel(N));

for i=1:numel(p)
    for j=1:numel(N)
        for k=1:numel(params)
            mesh = mkmesh_uniform([0, 1], N(j), p(i));
            master = mkmaster(mesh);

            % gauss points in physical coordinates
            dgg = master.phi'*mesh.dgnodes;

            [A,f] = assemble(mesh, params{k}{1}, params{k}{2}, params{k}{3}, params{k}{4});
            u = A \ f;
            ut = u(mesh.nn(:));
            un = reshape(ut, size(mesh.dgnodes));
            eg = (master.phi'*un - exact_solution{k}(dgg)); 
            J = master.dphi'*mesh.dgnodes;
            err2 = master.gw'*(J.*(eg.*eg));
            x = mesh.dgnodes(:);
            y = exact_solution{k}(x);
             figure;
             hold all;
             plot(x, ut, '.b', 'DisplayName', 'DG');
             plot( x, y, 'DisplayName', 'Exact');
             hold off;
             legend(gca, 'show', 'location', 'NorthWest');
            l2err = sqrt(sum(err2));
            errors(k, i, j) = l2err;
        end
    end
end

for k=1:numel(params)
    h = figure;
    set(h, 'units', 'inches', 'position', [1 1 4 4])
    set(h, 'PaperUnits','centimeters');
    set(h, 'Units','centimeters');
    pos=get(h,'Position');
    set(h, 'PaperSize', [pos(3) pos(4)]);
    set(h, 'PaperPositionMode', 'manual');
    set(h, 'PaperPosition',[0 0 pos(3) pos(4)]);

    hold all;
    for i=1:numel(p)
        y = squeeze(errors(k, i, :));
        loglog(N, y, 'DisplayName', sprintf('%s - Order %d', params{k}{5}, p(i)));
        pf = polyfit(log(N), log(y), 1);
        fprintf('%s - Order %d has rate %g.\n', params{k}{5}, p(i), pf(1));
    end
    hold off;
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    axis square;
    grid on;
    legend(gca, 'show', 'location', 'NorthEast');

    print(sprintf('../report/convergence_%s.pdf', params{k}{5}), '-dpdf');
end
