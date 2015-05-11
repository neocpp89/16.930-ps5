clear all; close all;

mesh = mkmesh_uniform([0, 1], 10, 2);
master = mkmaster(mesh);

[A,f] = assemble(mesh, 1, 0, 0, @(x) x.*x)
u = A \ f;
ut = u(1:end-2);
x = linspace(0,1,numel(ut));
y  = (-x.^3 + 7*x) / 6;
plot(x, ut, x, y)
