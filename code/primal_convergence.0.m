clear all; close all;

plot_solutions = 1;
% Ndof = 2.^[4:7]';
Ndof = 2.^[2:8]';
p = [1, 2, 3, 4, 5]';
Ndof_actual = zeros(numel(p), numel(Ndof));

params = { {1, 0, 0, @(x) x.*x, 'Poisson'}, ...
    {1e-2, 0, 1, @(x) 0*x, 'Convection-Diffusion'}, ...
    {1e-4, 1, 0, @(x) 0*x, 'Reaction-Diffusion'} };
% params = { {1e-2, 1, 0, @(x) -20*1e-2*(x.^3) + 1*(x.^5), 'Reaction-Diffusion'} };
% params = { {1e-2, 1, 0, @(x) 0*x, 'Reaction-Diffusion'} };
exact_solution = { @(x) x.*(13 - x.^3)/12, ...
   @(x) (1-exp(100*x)) ./ (1 - exp(100)), ...
   @(x) (exp(100*x)-exp(-100*x)) ./ (exp(100) - exp(-100)) };
% exact_solution = { @(x) x.^5 };
% exact_solution = { @(x) (exp(10*x)-exp(-10*x)) ./ (exp(10) - exp(-10)) };
errors = zeros(numel(params), numel(p), numel(Ndof));

for i=1:numel(p)
    for j=1:numel(Ndof)
        for k=1:numel(params)
            N = ceil(Ndof(j) / (1 + p(i)));
            mesh = mkmesh_uniform([0, 1], N, p(i));
            Ndof_actual(i, j) = numel(mesh.dgnodes);
            master = mkmaster(mesh);

            % gauss points in physical coordinates
            dgg = master.phi'*mesh.dgnodes;

            [A,f] = assemble(mesh, params{k}{1}, params{k}{2}, params{k}{3}, params{k}{4});
            u = A \ f;
            ut = u(mesh.nn(:));
            un = reshape(ut, size(mesh.dgnodes));
            eg2 = (master.phi'*un - exact_solution{k}(dgg)).^2; 
            J = master.dphi'*mesh.dgnodes;
            err2 = master.gw'*(J.*(eg2));
            x = mesh.dgnodes(:);
            y = exact_solution{k}(x);
            if (plot_solutions == 1)
                figure;
                hold all;
                plot(x, ut, '.b', 'DisplayName', 'DG');
                plot( x, y, 'DisplayName', 'Exact');
                title(sprintf('%s - Order %d (N = %d)', params{k}{5}, p(i), Ndof_actual(j)));
                hold off;
                legend(gca, 'show', 'location', 'Best');
            end
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
        loglog(Ndof_actual(i, :)', y, 'DisplayName', sprintf('%s - Order %d', params{k}{5}, p(i)));
        % pf = polyfit(log(N), log(y), 1);
        % only pick the last couple of points to do convergence rate
        pf = polyfit(log(Ndof_actual(i, end-1:end))', log(y(end-1:end)), 1);
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
