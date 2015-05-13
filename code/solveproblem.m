clear all; close all;

mesh = mkmesh_uniform([0, 1], 10, 1);
master = mkmaster(mesh);
mesh.bcv(1) = 0;

% [A,f] = assemble(mesh, 1, 0, 0, @(x) x.*x);
% [A,f] = assemble(mesh, 1e-1, 0, 0.5, @(x) x);
nu = 1;
b = 0;
c = 0;
% fn = @(x) -2*nu + b*x.*x;
fn = @(x) x.*x;
[A,f] = assemble(mesh, nu, b, c, fn);
u = A \ f;
ut = u(mesh.nn(:));
x = mesh.dgnodes(:);
y = x.*(13 - x.^3)/12;
% y = x.*x;
% yfn_reaction = @(r, x) (exp(r*x)-exp(-r*x)) ./ (exp(r) - exp(-r));
% y = yfn_reaction(sqrt(b/nu), x);
% yfn = @(x) x.*x;
% y = yfn(x);
plot(x, ut, 'xr-', x, y)