function [master] = mkmaster(mesh, porder)
if nargin == 1
    % try to make sure we can integrate the mass matrix at least...
    master.porder = 2*mesh.porder;
else
    master.porder = porder;
end
ng = ceil((master.porder + 1) / 2);
master.plocal = mesh.plocal;
np = numel(master.plocal);
master.edge = [1; mesh.porder+1]; % left and right
[x, w] = golub_welsch(ng);
master.gp = (x+1)/2;
master.gw = w/2;
master.phi = zeros(np, ng);
master.dphi = zeros(np, ng);
master.dleft = zeros(np, 1);
master.dright = zeros(np, 1);
P = lagrange(master.plocal);
for i=1:np
    master.phi(i, :) = polyval(P(i, :), master.gp);
    master.dphi(i, :) = polyval(polyder(P(i, :)), master.gp);
    master.dleft(i) = polyval(polyder(P(i, :)), master.plocal(1));
    master.dright(i) = polyval(polyder(P(i, :)), master.plocal(end));
end
master.P = P;
master.mass = master.phi * diag(master.gw) * master.phi';
wf = zeros(np, 1);
wf(1) = -0.5;
master.al = master.mass \ wf; % coefficients for face lifting operator rf(1) (left)

wf = zeros(np, 1);
wf(end) = -0.5;
master.ar = master.mass \ wf; % coefficients for face lifting operator rf(1) (right)
end
