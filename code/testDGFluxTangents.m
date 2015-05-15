classdef testDGFluxTangents < matlab.unittest.TestCase
    % Tests tangent matrices for convection-diffusion-reaction problems 

    methods (Test)
        function testFaceDGFluxMasterP1(testCase)
            % single element of order 1
            mesh = mkmesh_uniform([0, 1], 1, 1);
            master = mkmaster(mesh);

            jacobians = master.dphi'*mesh.dgnodes;

            % should be the same as the master
            testCase.verifyEqual(jacobians, ones(size(jacobians)), 'Abstol', 1e-10);

            c = 1;
            A = face_dg_flux_tangent(master, jacobians, c);
            testCase.verifyEqual(A, 0.5*[c+abs(c), c-abs(c); -c-abs(c), -c+abs(c)], 'Abstol', 1e-10);

            c = -1;
            A = face_dg_flux_tangent(master, jacobians, c);
            testCase.verifyEqual(A, 0.5*[c+abs(c), c-abs(c); -c-abs(c), -c+abs(c)], 'Abstol', 1e-10);
        end

        function testFaceDGFluxP1(testCase)
            % single element of order 1
            mesh = mkmesh_uniform([0, 0.5], 1, 1);
            master = mkmaster(mesh);

            jacobians = master.dphi'*mesh.dgnodes;

            % should be the half the size of master
            testCase.verifyEqual(jacobians, 0.5*ones(size(jacobians)), 'Abstol', 1e-10);

            % but dg inviscid flux tangents do not change
            c = 1;
            A = face_dg_flux_tangent(master, jacobians, c);
            testCase.verifyEqual(A, 0.5*[c+abs(c), c-abs(c); -c-abs(c), -c+abs(c)], 'Abstol', 1e-10);

            c = -1;
            A = face_dg_flux_tangent(master, jacobians, c);
            testCase.verifyEqual(A, 0.5*[c+abs(c), c-abs(c); -c-abs(c), -c+abs(c)], 'Abstol', 1e-10);
        end
    end
end
