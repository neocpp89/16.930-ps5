function [master] = mkmaster(mesh)
master.order = mesh.porder;
ng = ceil((master.order + 1) / 2);
master.plocal = mesh.plocal;
master.edge = [1; mesh.porder+1];
[x, w] = golub_welsch(ng);
master.gp = (x+1)/2;
master.gw = w/2;
% master.phi = ;
% master.dphi = ;
end
