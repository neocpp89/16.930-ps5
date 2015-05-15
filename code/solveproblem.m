clear all; close all;

mesh = mkmesh_uniform([0, 1], 100, 1);
master = mkmaster(mesh);
mesh.bcv(1) = 0;

% [A,f] = assemble(mesh, 1, 0, 0, @(x) x.*x);
% [A,f] = assemble(mesh, 1e-1, 0, 0.5, @(x) x);
nu = 1e-4;
b = 1;
c = 0;
fn = @(x) -2*nu + b*x.*x;
% fn = @(x) 0*x;
[A,f] = assemble(mesh, nu, b, c, fn);
u = A \ f;
ut = u(mesh.nn(:));
x = mesh.dgnodes(:);
% y = x.*(13 - x.^3)/12;
y = x.*x;
yfn_reaction = @(r, x) (exp(r*x)-exp(-r*x)) ./ (exp(r) - exp(-r));
yfn_convection = @(r, x) (exp(r*x)-1) ./ (exp(r) - 1);
% y = yfn_convection(c/nu, x);
xx = linspace(0,1,200);
yy = yfn_reaction(sqrt(b/nu), xx);
yfn = @(x) x.*x;
% y = yfn(x);
yy= yfn(xx);
plot(x, ut, 'xr-', xx, yy)