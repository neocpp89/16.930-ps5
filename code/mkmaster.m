function [master] = mkmaster(mesh)
master.order = mesh.porder;
master.plocal = mesh.plocal;
master.edge = [1; mesh.porder];
master.gw = 1;
master.gp = 1; 
end
