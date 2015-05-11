clear all; close all;

mesh = mkmesh_uniform([0, 1], 1000, 3);
master = mkmaster(mesh);

[A,f] = assemble(mesh, 1, 0, 0, @(x) x.*x);
% [A,f] = assemble(mesh, 0, 0, 0.5, @(x) x);
u = A \ f;
ut = u(mesh.nn(:));
x = mesh.dgnodes(:);
y = x.*(13 - x.^3)/12;
% y = x.*x;
plot(x, ut, '.b', x, y)
